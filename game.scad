// Game organizer library

include <cards.scad>

// TODO: hex boxes and grids from calico-box and civ-box

// naming conventions
// A  angle
// C  color
// D  depth, diameter, thickness
// H  height
// N  number
// P  polygon (list of points)
// Q  quality setting (draft, final)
// R  radius
// S  spin (rotation vector)
// V  vector  [W, H] or [W, D, H] or [x, y, z]
// W  width

Qdraft = 15;  // 24 segments per circle (aligns with axes)
Qfinal = 5;
$fa = Qdraft;
$fs = 0.1;

EPSILON = 0.01;
MICRON = 0.001;
PHI = (1+sqrt(5))/2;

// printer & filament metrics
Vprinter = [180, 180, 180];
Hflayer = 0.25;
Dfwidth = 0.70;  // extrusion width
Dfoverlap = Hflayer * (1 - PI/4);  // overlap between paths
Dfpath = tfloor(Dfwidth - Dfoverlap);  // width multiplier for walls
echo(Vprinter=Vprinter, Hflayer=Hflayer, Dfwidth=Dfwidth,
     Dfoverlap=mround(Dfoverlap), Dfpath=Dfpath);

// organizer metrics
Dwall = 2.0;
Hfloor = Dwall;
Dthick = 3.0;  // for heavier, stiffer walls
Dthin = 1.0;  // for thin divider walls
Dgap = Dfpath/4;
Dcut = eround(Dwall/3);  // cutting margin for negative spaces
Djoiner = EPSILON;  // overlap margin for joining parts
echo(Dwall=Dwall, Hfloor=Hfloor, Dgap=Dgap, Dcut=Dcut, Djoiner=Djoiner);
Rext = 3.0;  // external corner radius
Rint = Rext - Dwall;  // internal corner radius
echo(Rext=Rext, Rint=Rint);
Hstack = 1.0;  // stacker joint height
Drim = 1.2;  // stacker rim thickness
echo(Hstack=Hstack, Drim=Drim);
Avee = 60;  // default angle for notches (TODO: replace with Atab)
Atab = 60;  // default angle for tabs & notches
Ahex = 60;  // default angle for hexagons & triangles
Adeep = 75;  // default angle for deep notches
Arack = 75;  // default angle for card & tile racks
Adraw = 3;  // default slope for draw trays (TODO: are these obsolete?)
Sup = [90, 0, 0];
Sdown = [-90, 0, 0];
echo(Avee=Avee, Ahex=Ahex, Atab=Atab, Arack=Arack, Adraw=Adraw);
Dthumb = 25.0;  // index hole diameter
echo(Dthumb=Dthumb);

// game box interior
Vgame = [288, 288, 69];  // typical FFG box interior
Hwrap = 55;  // cover art wrap ends here, approximately
echo(Vgame=Vgame, Hwrap=Hwrap);

// component metrics
Nplayers = undef;  // number of players (for mats & other per-player items)
Hboard = 2.5;  // thickness of cardboard & similar flat components
Hmat = Hboard;  // mats: trackers, player boards, holding areas
Htile = Hboard;  // tiles: hexes, maps, plaques
Htoken = Hboard;  // tokens: coins, points, units
// hex tiles
Rhex = Dthumb;  // hex major radius = side length = grid spacing
Rhex_grid = Rhex;  // size of gridded hexes (may overflow spacing)
Rhex_single = Rhex;  // size of ungridded hex tiles
echo(Rhex=Rhex, Rhex_grid=Rhex_grid, Rhex_single=Rhex_single);
// chips & counters
Dchip = 40.0;
Rchip = Dchip / 2;
Hchip = 3.4;
// dice
Ddice = 16;

// available space
Hmanual = 1.0;
Hceiling = Vgame.z - eceil(Hmanual, 0.5);  // vertical space under manuals
Hmain = Hceiling;  // vertical space under manuals and boards
echo(Hmanual=Hmanual, Hceiling=Hceiling, Hmain=Hmain);

// container metrics
Rfoot = Rint - Dgap;  // concentric with Rint & Rext with nesting gap
Hfoot = 1.0;  // TODO: replace with Hbase or Hstack?
Hhead = 1.0;  // headroom above contents
Hlip = Hhead + Hfoot;  // wall height above contents, scoops, etc.
Htab = 1 - Hflayer;
Htray = 13.0;
Vtray = [72, 100, Htray];
Vfoot = volume(Vtray/8, Hfoot);
echo(Vtray=Vtray, Htray=Htray, Hlip=Hlip);
echo(Vfoot=Vfoot, Hfoot=Hfoot, Rfoot=Rfoot);
Vbox = Vtray;  // TODO
echo(Vbox=Vbox);

// minimum sizes and rounding
function eround(x, e=EPSILON) = e * round(x/e);
function eceil(x, e=EPSILON) = e * ceil(x/e);
function efloor(x, e=EPSILON) = e * floor(x/e);
function tround(x) = eround(x, e=0.05);  // twentieths of a millimeter
function tceil(x) = eceil(x, e=0.05);  // twentieths of a millimeter
function tfloor(x) = efloor(x, e=0.05);  // twentieths of a millimeter
function mround(x) = eround(x, e=MICRON);  // microns
function mceil(x) = eceil(x, e=MICRON);  // microns
function mfloor(x) = efloor(x, e=MICRON);  // microns
function lround(x) = eround(x, e=Hflayer);  // layers
function lceil(x) = eceil(x, e=Hflayer);  // layers
function lfloor(x) = efloor(x, e=Hflayer);  // layers
function pround(x) = eround(x, e=Dfpath);  // paths
function pceil(x) = eceil(x, e=Dfpath);  // paths
function pfloor(x) = efloor(x, e=Dfpath);  // paths

// tidy measurements
function vround(v) = [for (x=v) tround(x)];
function vceil(v) = [for (x=v) tceil(x)];
function vfloor(v) = [for (x=v) tfloor(x)];

// max & min for geometric quantities: keep sign but order by magnitude
function absmax(x, y) = let (ax=abs(x), ay=abs(y))
    ax < ay ? y : ay < ax ? x : max(x, y);
function absmin(x, y) = let (ax=abs(x), ay=abs(y))
    ax < ay ? x : ay < ax ? y : min(x, y);
// normalized area & volume vectors
function area(size, wide=undef) =
let (v=is_list(size) ? [size.x, size.y] : [size, size]) [
    // calculate area with optional wide or tall override
    is_undef(wide) ? v.x : wide ? absmax(v.x, v.y) : absmin(v.x, v.y),
    is_undef(wide) ? v.y : wide ? absmin(v.x, v.y) : absmax(v.x, v.y),
];
function volume(size, height=undef, wide=undef) =
let (v=is_list(size) ? [size.x, size.y, size.z] : [size, size, size]) [
    // calculate volume with optional height, wide, or tall override
    is_undef(wide) ? v.x : wide ? absmax(v.x, v.y) : absmin(v.x, v.y),
    is_undef(wide) ? v.y : wide ? absmin(v.x, v.y) : absmax(v.x, v.y),
    is_undef(height) ? v.z : height,
];
function grid(size=1) =
let (v=is_list(size) ? [size.x, size.y] : size ? [1, size] : [0, 0]) [
    // scalar: [1, depth]; 1-vector: [width, 1]; 2-vector: [width, depth]
    is_undef(v.x) ? 1 : v.x,
    is_undef(v.y) ? 1 : v.y,
];

// deck & box dimensions
function deck_volume(n=1, size=Vcard, height=Hcard) =
    volume(size, n * height);
// box volume is equal to deck volume plus:
//     2*Rext in the X & Y dimensions
//     Hfloor+lip in the Z dimenions
// width overrides the X dimension if set (with no added margin)
function deck_box_volume(n=0, size=Vcard, height=Hcard, width=0, lip=Hlip) =
    let (v = deck_volume(n, size, height),
         w = width ? width : v.z ? v.z + 2*Rext : Vtray.x,
         d = v.y + 2*Rext,
         h = v.x + Hfloor + lip)
        // echo(v=v, w=w, d=d, h=h)
        // echo(cards=floor((w-2*Rext)/height))
        // echo(piles=floor((w-2*Rext)/height)/12)
        [w, d, h];

// utility functions
function sum(v) = v ? [for(p=v) 1]*v : 0;
function swapxy(v) = [v.y, v.x, if (2<len(v)) for (i=[2:len(v)-1]) v[i]];
function unit_axis(n) = [for (i=[0:2]) i==n ? 1 : 0];
function numeric_flag(x, default) = is_num(x) ? x : x ? default : 0;

// transformations
module colorize(c=undef, alpha=undef) {
    // skip the color() call if both parameters are undef
    if (is_undef(c) && is_undef(alpha)) children();
    else color(c, alpha) children();
}
module flatten(size, height=undef, space=undef, angle=undef) {
    // shear and flatten with fixed sides (like flattening a cardboard box)
    v = volume(size, height);
    c = is_undef(space) ? v.x : max(space, v.x);
    dx = min(c - v.x, v.z - EPSILON);
    A = is_undef(angle) ? acos(dx/v.z) : max(angle, EPSILON);
    x = v.x;
    z = v.z*sin(A);
    xc = x + v.z*cos(A);
    mlean = [
        [1, 0, -cos(A), xc-x],
        [0, 1, 0, 0],
        [0, 0, sin(A), 0],
    ];
    multmatrix(m=mlean) children();
}
module lean(size, height=undef, space=undef, angle=undef) {
    // shear and rotate with fixed volume (like leaning cards against a box)
    v = volume(size, height);
    c = is_undef(space) ? v.x : max(space, v.x);
    function solve() = let (
        // x^4 + Bcx^3 + Cx^2 + E = 0, via Ferrari's method
        B = -2*c,
        C = c^2 - v.z^2,
        E = v.x^2*v.z^2,
        a = -3*B^2/8 + C,
        b = B^3/8 - B*C/2,
        g = -3*B^4/256 + C*B^2/16 + E,
        p = -a^2/12 - g,
        q = -a^3/108 + a*g/3 - b^2/8,
        r = -q/2 + sqrt(q^2/4 + p^3/27),
        u = r^(1/3),
        y = -5/6*a + u - p/(3*u),
        w = sqrt(a + 2*y),
        v = sqrt(-(3*a + 2*y + 2*b/w)),
        x = -B/4 + (w-v)/2)
        x;
    x = is_undef(angle) ? solve() : v.x/sin(max(angle, EPSILON));
    A = is_undef(angle) ? asin(v.x/x) : max(angle, EPSILON);
    z = v.z*sin(A);
    xc = x + v.z*cos(A);
    mshear = [
        [1, 0, 0, -v.x/2],
        [0, 1, 0, 0],
        [-1/tan(A), 0, 1, v.x/tan(A)/2],
    ];
    mrotate = [
        [sin(A), 0, -cos(A), xc-x/2],
        [0, 1, 0, 0],
        [cos(A), 0, sin(A), 0],
    ];
    multmatrix(m=mrotate) multmatrix(m=mshear) children();
}
module raise(z=Hfloor+EPSILON) {
    translate([0, 0, z]) children();
}

// 3D shapes
module fillet(rint=undef, rext=undef) {
    // round inside corners to radius rint, outside corners to rext
    // parameter order reflects typical usage of "fillet" for inside rounding
    if (rint && rext)
        fillet(rext=rext) fillet(rint=rint) children();
    else if (rint)
        offset(r=-rint) offset(delta=rint) children();
    else if (rext)
        offset(r=rext) offset(delta=-rext) children();
    else children();
}
module rounded_square(size, r=Rext) {
    // creates a rounded square with corners of radius r
    v = area(size);
    if (min(v)/2 <= r) stadium_fill(v);
    else fillet(rext=r) square(v, center=true);
}
module quarter_round(size, r=Rext) {
    // creates a square in the first quadrant with rounded outside corner
    v = area(size);
    o = v - [r, r];
    if (r) hull() {
        square([o.x, v.y]);
        square([v.x, o.y]);
        translate(o) intersection() {
            circle(r=r);
            square(r);
        }
    } else square(size);
}
module stadium(h, r=undef, d=undef) {
    // creates a stadium with rectangle height h and radius r,
    // centered on the Y axis
    radius = abs(is_undef(d) ? r : d/2);
    height = abs(h);
    hull() {
        if (height) square([2*radius, height], center=true);
        for (i=[-1,+1]) translate([0, i*height/2]) circle(radius);
    }
}
module stadium_fill(size) {
    // creates a stadium sized to fit the given area
    v = area(size, wide=false);  // aligned vertically for stadium module
    v0 = area(size);  // original alignment
    r = abs(v.x)/2;
    h = abs(v.y) - 2*r;
    a = (v.x == v0.x ? 0 : -90);
    rotate(a) stadium(h, r);
}
module semistadium(h, r=undef, d=undef) {
    // creates a semistadium with rectangle height h and radius r,
    // centered on the positive Y axis
    radius = abs(is_undef(d) ? r : d/2);
    s = h ? sign(h) : 1;  // there's no negative 0, so default to positive
    hull() {
        if (h) translate([0, h/2]) square([2*radius, abs(h)], center=true);
        translate([0, h]) intersection() {
            circle(radius);
            translate([0, s*radius]) square(2*radius, center=true);
        }
    }
}
module semistadium_fill(size) {
    // creates a semistadium sized to fit the given area,
    // centered on the positive Y axis, scaled to width if necessary
    v = area(size);
    r = min(abs(v.x/2), abs(v.y));
    h = max(abs(v.y) - r, 0);
    s = h ? sign(v.y) : [v.x / r / 2, sign(v.y)];
    scale(s) semistadium(h, r);
}
module capsule(h, r=undef, d=undef) {
    // creates a capsule with cylinder height h and radius r,
    // centered on the Z axis
    radius = abs(is_undef(d) ? r : d/2);
    height = abs(h);
    hull() {
        if (height) cylinder(h=h, r=radius, center=true);
        for (i=[-1,+1]) translate([0, 0, i*height/2]) sphere(radius);
    }
}
module semicapsule(h, r=undef, d=undef) {
    // creates a semicapsule with cylinder height h and radius r,
    // centered on the positive Z axis
    radius = abs(is_undef(d) ? r : d/2);
    s = h ? sign(h) : 1;  // there's no negative 0, so default to positive
    hull() {
        if (h) translate([0, 0, h/2])
            cylinder(h=abs(h), r=radius, center=true);
        translate([0, 0, h]) intersection() {
            sphere(radius);
            scale([1, 1, s]) cylinder(h=2*radius, r=2*radius);
        }
    }
}

module prism(size=undef, height=undef, r=undef, rint=undef, rext=undef,
             scale=1, center=false) {
    v = is_undef(size) ? undef : volume(size, height);
    h = is_undef(height) ? v.z : height;
    ri = is_undef(rint) ? is_undef(r) ? 0 : r : rint;  // inside turns
    re = is_undef(rext) ? is_undef(r) ? 0 : r : rext;  // outside turns
    linear_extrude(height=h, scale=scale, center=center) fillet(ri, re) {
        if (is_undef(v)) children();
        else square(area(v), center=true);
    }
}
module box_frame(size=Vgame, height=undef, wall=Dwall, wrap=Hwrap, gap=Dgap) {
    // create the outline of a box with given interior and thickness
    vint = volume(size, height);
    vext = vint + [2*wall, 2*wall, wall];
    dwall = wall - gap;  // shrink the wall to leave a small gap
    vcut = vext - [2*dwall, 2*dwall, 2*dwall];
    raise(vint.z - vext.z) difference() {
        prism(vext);
        raise(dwall) {
            for (n=[0:2]) for (i=[-1,+1])  // sides
                translate(i*unit_axis(n) * 4/3*dwall) prism(vcut);
        }
    }
    if (is_num(wrap)) {
        raise(wrap) linear_extrude(height=2*wall) difference() {
            square([vext.x, vext.y], center=true);
            square([vcut.x, vcut.y], center=true);
        }
    }
}

module floor_thumb_cut(size, height=undef, d=Dthumb, r=Rext, mirror=false, cut=Dcut) {
    v = volume(size, height);
    dy = d/2;  // depth of thumb round
    s = mirror ? [-1, +1] : [+1];
    h = v.z + 2*cut;
    raise(-cut) {
        // thumb round
        for (s=s) scale([1, s]) translate([0, -cut-v.y/2])
            prism(height=h, rint=r) {
                // approximate width of opening at the tangents
                // (quantization of $fa causes a small difference from ideal)
                axis = d/2 + r;
                span = 2*axis*cos(asin(r/axis));
                translate([0, cut/2]) square([span, cut-EPSILON], center=true);
                semistadium(dy - d/2 + cut, d=d);
            }
        // bottom index hole
        prism(height=h) circle(d=d);

    }
}
module wall_vee_cut(size, height=undef, angle=Avee, cut=Dcut, fillet=true) {
    a0 = max(EPSILON, min(angle, 90));
    v = volume(size, height);
    run = a0 < 90 ? 1/tan(a0) : 0;
    a1 = 90 - a0/2;
    y1 = v.z - Rext;
    y2 = v.z;
    y3 = v.z + cut;
    x0 = v.x/2;
    x1 = x0 + y1*run;
    x2 = x0 + y2*run;
    x3 = x2 + Rext/tan(a1);
    d = v.y + 2*cut;
    if (fillet) {
        ptop = [
            [x3, y3], [x3, y2], [x2, y2], [x1, y1],
            [-x1, y1], [-x2, y2], [-x3, y2], [-x3, y3],
        ];
        pbot = [[x2, y2], [x0, 0], [-x0, 0], [-x2, y2]];
        rotate([90, 0, 0]) {
            prism(height=d, rint=Rext, center=true) polygon(ptop);
            prism(height=d, rext=Rint, center=true) polygon(pbot);
        }
    } else {
        pcut = [[x2, y3], [x2, y2], [x0, 0], [-x0, 0], [-x2, y2], [-x2, y3]];
        rotate([90, 0, 0]) prism(height=d, center=true) polygon(pcut);
    }
}
module hex_cut(size, height=undef, cut=Dcut) {
    wall_vee_cut(size=size, height=height, angle=Ahex, cut=cut, fillet=false);
}

module deck_box(n=0, size=Vcard, height=Hcard, width=0, lip=Hlip, draw=false,
                feet=false, color=undef) {
    vbox = deck_box_volume(n=n, size=size, height=height, width=width);
    shell = area(vbox);
    well = shell - 2 * area(Dwall);
    hole = shell - 2 * area(shell.y/5);
    echo(vbox=vbox);
    translate([vbox.x/2, 0]) colorize(color) difference() {
        // outer shell
        prism(vbox, r=Rext);
        // card well
        raise(Hfloor) prism(well, height=vbox.z, r=Rint);
        // base round (if it fits)
        dh = min(hole.x, hole.y);  // hole diameter
        if (3/5*Dthumb <= dh && dh <= Dthumb)
            raise(-Dgap) prism(height=vbox.z) stadium_fill(hole);
        else raise(-Dgap) prism(hole, height=vbox.z, r=Dthumb/2);
        dtop = vbox.y - 4*Rext;  // maximum notch width
        if (draw) {
            // thumb cut
            vthumb = [Dthumb/sin(Avee), 2*Dwall, Dthumb];
            translate([(Dwall-vbox.x)/2, 0, vbox.z-vthumb.z])
                rotate(90) wall_vee_cut(vthumb);
            // front cut
            adraw = Adeep;
            hvee = vbox.z - Hfloor;  // maximum height
            dxvee = hvee / tan(adraw);
            vdraw = [dtop - 2*dxvee, 2*Dwall, hvee];
            translate([(vbox.x-Dwall)/2, 0, Hfloor])
                rotate(90) wall_vee_cut(vdraw, angle=adraw);
        } else {
            // side cuts
            zvee = min(vbox.z/2, dtop*sin(Avee)/2);
            hvee = vbox.z-zvee;
            xvee = tround(zvee/sin(Avee));
            vend = [xvee, vbox.x, zvee];
            echo(vbox=vbox, dtop=dtop, zvee=zvee, xvee=xvee);
            raise(hvee) rotate(90) wall_vee_cut(vend);  // end vee
        }
    }
    // feet
    if (feet) colorize(color) for (i=[-1,+1]) {
        // center feet in the available space
        yin = Dthumb/sin(Avee) + Rext/tan(Avee);
        yout = vbox.y/2 - Rext;
        yfoot = (yin + yout) / 2;
        translate([Rext-Dwall, i*yfoot, vbox.z-Rext-Rint]) intersection() {
            translate([-3/2*Rext, 0]) cube(3*Rext, center=true);
            sphere(Rext);
        }
    }
    translate([vbox.x + Dgap, 0]) children();
}
module draw_box(n=0, size=Vcard, height=Hcard, width=0, lip=Hlip, feet=true,
                color=undef) {
    deck_box(n=n, size=size, height=height, width=width, lip=lip, draw=true,
             feet=feet, color=color) children();
}

module tray_feet_cut(size=Vtray, height=undef, foot=Vfoot) {
    if (foot.z) {
        v = volume(size, height);
        d = Rext - Rfoot;  // margin between foot and tray
        o = (area(v) - area(foot))/2 - area(d);
        for (i=[-1,+1]) for (j=[-1,+1])
            translate([i * o.x, j * o.y])
                tray_foot(cut=Dcut);
    }
}
module tray_foot(size=Vfoot, height=undef, r=Rfoot, cut=0) {
    // creates feet for nesting trays, or set cut=Dcut to make the leg socket
    vfoot = volume(size, height);
    vslot = volume(vfoot, Hfloor/2) - volume(2*r, 0);
    vleg = vslot - volume(Dgap, Hflayer);  // fit tolerance + room for glue
    if (cut) {
        raise(-cut) prism(vslot, height=cut+vslot.z);
        %raise(-vfoot.z) tray_foot();
    } else {
        prism(vfoot, r=r);
        raise(vfoot.z-EPSILON) prism(vleg, height=vleg.z+EPSILON);
    }
}
module card_well(size=Vtray, height=undef, cut=Dcut) {
    vtray = volume(size, height);
    vwell = volume(area(vtray) - 2*area(Dwall), vtray.z-Hfloor);
    raise(Hfloor) {
        // card well
        prism(vwell, height=vwell.z+cut, r=Rint);
        // thumb vee
        span = Dthumb + 2*Rint;
        dmax = (vwell.x - span) / 4;  // maximum spread of vee at top
        amin = atan((vtray.z-Hfloor)/dmax);  // minimum vee angle
        echo(span=span, dmax=dmax, amin=amin);
        angle = max(Avee, eround(amin, Qfinal));
        translate([0, Dwall-vtray.y]/2)
            wall_vee_cut([span, Dwall, vtray.z-Hfloor], angle=angle, cut=cut);
    }
    floor_thumb_cut(vtray, cut=cut);
}
module card_tray(size=Vtray, height=undef, cards=0, feet=true, color=undef) {
    vtray = volume(size, height);
    colorize(color) difference() {
        prism(vtray, r=Rext);
        card_well(vtray);
        if (feet && Hfoot) tray_feet_cut(vtray);
    }
    %raise()  // card stack
        if (cards) deck(cards) children();
        else children();
}
module draw_tray(size=Vtray, height=undef, slope=Adraw, color=undef) {
    vtray = volume(size, height);
    vwell = area(vtray) - 2*area(Dwall);
    hface = tan(slope) * vwell.y + Hfloor;
    mslope = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, -sin(slope), 1, 0],
    ];
    colorize(color) difference() {
        prism(vtray, r=Rext);
        floor_thumb_cut(vtray);
        tray_feet_cut(vtray);
        // sloped floor
        raise(Hfloor/2 + hface/2) multmatrix(m=mslope)
            prism(vwell, height=vtray.z+Dcut, r=Rint);
        // open front
        translate([0, -vtray.y/2, hface]) rotate(90)
            wall_vee_cut([2*Rext, vtray.x + Dcut, vtray.z-hface]);
    }
}
module deck(n=10, size=Vcard, height=Hcard, up=false, color=undef) {
    v = deck_volume(n=n, size=size, height=height);
    spin = up ? [0, 90, 0] : 0;
    lift = up ? v.x/2 : 0;
    raise(lift) rotate(spin) colorize(color) prism(v, height=v.z);
    translate(spin ? [v.z, 0, 0] : [0, 0, v.z]) children();
}
module creasing_tool(n=10, size=Vcard, height=Hcard, color=undef) {
    // block for creasing index wrappers
    v = deck_volume(n=n, size=size, height=height);
    echo(creasing_tool=v);
    colorize(color) prism(height=v.z) difference() {
        square(area(v), center=true);
        circle(d=Dthumb);
    }
}
module deck_divider(size=Vcard_divider, height=Hcard_divider,
                    up=false, color=undef) {
    // vertical divider for wide deck boxes
    v = volume(size, height);
    spin = up ? [0, 90, 0] : 0;
    lift = up ? v.x/2 : 0;
    raise(lift) rotate(spin) colorize(color) {
        xthumb = 2/3 * Dthumb;  // depth of indentation
        y0 = v.y/2;
        y1 = xthumb/sin(Avee);
        y2 = y1/2;
        x0 = v.x/2;
        x1 = x0 - xthumb;
        poly = [
            [+x0, -y0], [+x0, -y1], [+x1, -y2],
            [+x1, +y2], [+x0, +y1], [+x0, +y0],
            [-x0, +y0], [-x0, +y1], [-x1, +y2],
            [-x1, -y2], [-x0, -y1], [-x0, -y0],
        ];
        prism(height=v.z, r=Rext) polygon(poly);
    }
    translate(spin ? [v.z, 0, 0] : [0, 0, v.z]) children();
}
module tray_divider(size=Vcard_divider, height=Hcard_divider,
                    index=Vtray.y, color=undef) {
    // horizontal divider for card trays
    v = volume(size, height);
    colorize(color) prism(height=v.z, r=Rext) difference() {
        square(area(v), center=true);
        // match the index holes in the underlying card tray
        projection() floor_thumb_cut([v.x, index, v.z], r=0, mirror=true);
    }
    translate([0, 0, v.z]) children();
}
module grid_divider(size=Vtray, height=undef, grid=[2, 3], wall=Dthick,
                    color=undef) {
    // rectangular grid divider for boxes and trays
    v = volume(size, height);
    grid = area(grid);
    function section(n, x) = (x + wall) / n;
    cell = [section(grid.x, v.x), section(grid.y, v.y)];
    echo(grid=grid, v=v, cell=cell);
    %prism(v);
    colorize(color) {
        for (i=[1:grid.x-1]) translate([i*cell.x - v.x/2 - wall/2, 0])
            prism([wall, v.y, v.z]);
        for (i=[1:grid.y-1]) translate([0, i*cell.y - v.y/2 - wall/2])
            prism([v.x, wall, v.z]);
    }
}
module scoop_well(size, height=undef, rint=Rint, rscoop=2*Rext, lip=Hlip,
                  cut=Dcut) {
    v = volume(size, height);
    hmax = v.z - lip;  // leave room for nesting feet
    rmax = min(v.x, v.y) / 4;  // limit radiuses to safe values
    rn0 = min(rint, rmax);
    rn1 = min(rscoop, rmax);
    hull() {
        raise(v.z) prism(v, height=cut, r=rn0);
        for (angle=[0:$fa:90]) {
            cz = 1-cos(angle);
            cx = 1-sin(angle);
            htier = hmax * cz;
            vtier = v - 2 * cx * area(rn1);
            rmax = min(vtier)/2 - EPSILON;
            rtier = min(rmax, cx * (rn1 - rn0) + rn0);
            raise(htier) prism(vtier, height=v.z-htier+EPSILON, r=rtier);
        }
    }
}
module scoop_tray(size=Vtray, height=undef, grid=1, rscoop=2*Rext, lip=Hlip,
                  color=undef) {
    // tray with scoop bottom and optional grid of wells
    grid = area(grid);
    v = volume(size, height);
    function section(n, x) = (x - Dwall) / n;
    cell = [section(grid.x, v.x), section(grid.y, v.y)];
    colorize(color) difference() {
        prism(v, r=Rext);
        raise(Hfloor) for (i=[1/2:grid.x]) for (j=[1/2:grid.y]) {
            translate(area(v)/2 - area(Dwall)/2 - [i*cell.x, j*cell.y])
            scoop_well(cell - area(Dwall), height=v.z-Hfloor, rscoop=rscoop,
                       lip=lip);
        }
    }
}

module hex(points=[[0, 0]], rhex=undef, grid=Rhex_grid, merge=Djoiner) {
    rhex = is_undef(rhex) ? len(points) == 1 ? Rhex_single : grid : rhex;
    x1 = sin(Ahex) * rhex;
    y1 = rhex / 2;
    phex = [[0, rhex], [-x1, y1], [-x1, -y1], [0, -rhex], [x1, -y1], [x1, y1]];
    dx = sin(Ahex) * grid;
    dy = grid;
    offset(delta=-merge) offset(delta=merge) for (p=points) {
        translate([2 * (p.x + p.y/2) * dx, 1.5 * p.y * dy]) {
            polygon(phex);
        }
    }
}
module hex_tile(points=[[0, 0]], n=0, height=Htile, rhex=undef, grid=Rhex_grid) {
    h = n ? eceil(n * height) : height;
    prism(height=h) hex(points, rhex=rhex, grid=grid);
}
module hex_tray(points=[[0, 0]], n=0, height=Htile, rhex=undef, grid=Rhex_grid,
                lip=Hlip, hole=Dthumb) {
    h = n ? eceil(n * height) + Hfloor + lip : height;
    difference() {
        prism(height=h, rint=Rint, rext=Rext)
            offset(delta=Rext) hex(points, rhex=rhex, grid=grid);
        raise(Hfloor) prism(height=h, rint=Rext, rext=Rint)
            offset(delta=Rint) hex(points, rhex=rhex, grid=grid);
        if (hole)
            raise(-Dcut) prism(height=Hfloor+2*Dcut, r=Rext)
            hex(points, rhex=hole/2, grid=grid);
    }
}
module hex_box(points=[[0, 0]], n=0, height=Htile, rhex=undef, grid=Rhex_grid,
               depth=undef, r=Rext, wall=undef, hole=undef, thick=false,
               stack=false, base=undef, lip=undef, color=undef) {
    rhex = is_undef(rhex) ? len(points) == 1 ? Rhex_single : grid : rhex;
    hole = numeric_flag(hole, Dthumb);
    stack = numeric_flag(stack, Hstack);
    h = n ? lceil(n * height) + Hfloor + Hhead + stack : height;
    wall = is_undef(wall) ? wall_thickness(wall, thick) : wall;
    difference() {
        stacker(height=h, depth=depth, r=r, wall=wall, stack=stack, base=base, lip=lip)
            offset(delta=r) hex(points, rhex=rhex, grid=grid);
        if (hole) punch(Hfloor) fillet(rext=r-wall) hex(points, hole/2, grid);
    }
}

module chip_tray(n=20, rows=5, color=undef) {
    r = Rchip + Dgap;
    h = lround(5/6*r);  // depth of slot
    a = asin((r-h)/r);
    w = 2*r*cos(a);  // width of slot at surface
    overhang = r - w/2;
    rail = 2*overhang + Rint;
    slot = [w, Hchip*n + Rint, h];
    well = [(slot.x + rail) * rows - rail, slot.y, slot.z];
    v = well + [2*Dwall, 2*Dwall, Hfloor];
    colorize(color) difference() {
        prism(v, r=Rext);
        for (i=[0:rows-1]) {
            translate([slot.x/2 - well.x/2 + (slot.x + rail)*i, 0, r+Hfloor]) {
                rotate([90, 0, 0]) cylinder(r=r, h=slot.y, center=true);
                %rotate([90, 0, 0]) cylinder(r=Rchip, h=Hchip*n, center=true);
            }
        }
    }
}

module tile_rack(n, size, angle=Arack, margin=Rext, lip=Hlip, color=undef) {
    vtile = volume(size, wide=true);
    echo(vtile=vtile);
    width = n * size.x + 2*margin;  // total width
    // size (hypotenuse) of back and foot rests
    back = max(vtile.x/2, vtile.y) + margin;
    zback = round(back * sin(angle));
    yback = zback/tan(angle);
    height = zback + Hfloor;
    depth = lceil(yback + (vtile.z+Dgap)*sin(angle) + 2*margin);
    yfoot = depth - yback - 2*margin;
    zfoot = (yfoot)/tan(angle);
    foot = yfoot/sin(angle);
    zlip = lround(zfoot + lip);
    echo(back=back, foot=foot);
    echo(yback=yback, yfoot=yfoot);
    echo(zback=zback, zfoot=zfoot, zlip=zlip);
    echo(height=height, depth=depth);
    shell = [width, depth, height];
    colorize(color) difference() {
        prism(shell, r=margin);
        well = [width+2*Dcut, foot, back+Dcut];
        translate([-width/2-Dcut, margin-depth/2, zfoot+Hfloor]) hull() {
            cube(well);
            rotate([angle-90, 0, 0]) cube(well);
        }
        translate([-width/2, -depth/2-Dcut, zlip+Hfloor])
            cube([width+2*Dcut, margin+2*Dcut, height-zlip+Dcut]);
    }
    %raise(Hfloor/2) cube([width, depth, Hfloor], center=true);
    %for (n=[1:n])
        translate([n*vtile.x-width/2+margin, yfoot+margin-depth/2, Hfloor])
        rotate([90-angle, 0, 180]) cube([vtile.x, vtile.z, vtile.y]);
    %translate([vtile.y/2, yfoot+margin-depth/2, Hfloor])
        rotate([90-angle, 0, 180]) cube([vtile.y, vtile.z, vtile.x]);
}

// tabs & notches
function tab_angle(v, w, r) =
    echo(v=v, w=w, r=r)
    // find the widest angle that fits between the tab shoulders
    // https://math.stackexchange.com/a/4479659/88237
    let (dc = [v.x/2-w/2, v.y-r],  // widest shoulder position
         dt = sqrt(dc.x^2 + dc.y^2 - r^2))  // corner -> shoulder tangent
        atan((dc.x*r + dc.y*dt) / (dc.x*dt - dc.y*r));
module tab(size, w1=undef, w2=undef, angle=Atab, rint=Rint, rext=Rext,
           joiner=Djoiner) {
    // create a tab shape inside a given area
    //   size    maximum extent of tab, including base rounding
    //   w1      base width
    //   w2      top width
    //   angle   rise angle
    //   rint    base rounding
    //   rext    top rounding
    //   joiner  depth below baseline (for joining parts)
    v = area(size);
    // adjust the angle to fit the available space, if needed
    min_angle = w2 ? tab_angle(v, w2, rint) : EPSILON;
    angle = w1 && w2 ? atan2(v.y, (w1-w2)/2) : max(angle, min_angle);
    dx1 = rint/tan(90 - angle/2);  // distance x1-x0
    dx2 = v.y/tan(angle);  // distance x2-x1
    x1 = w1 ? w1/2 : w2 ? w2/2 + dx2 : v.x/2 - dx1;  // base corner
    x2 = w2 ? w2/2 : x1 - dx2;  // top corner
    x0 = x1 + dx1;  // base tangent
    xmax = max(x0, x2);  // widest point
    xt = xmax + rext;  // base turnaround
    y0 = -2*rext - EPSILON;  // bottom of turnaround
    p = [
        [x2, v.y], [x1, 0], [xt, 0], [xt, y0],
        [-xt, y0], [-xt, 0], [-x1, 0], [-x2, v.y],
    ];
    echo(a=angle, w0=mround(2*x0), w1=mround(2*x1), w2=mround(2*x2));
    intersection() {
        fillet(rint, rext) polygon(p);
        translate([0, v.y/2]) square([2*xmax, v.y+2*joiner], center=true);
        translate([0, v.y/2]) square([v.x, v.y+2*joiner], center=true);
    }
}
module hex_tab(size=undef, rhex=undef, angle=Ahex, r=Rext, joiner=Djoiner) {
    ws = r/tan(90 - angle/2);  // shoulder width
    v = area(is_undef(size) ? 2*rhex + 2*ws : size);  // safe default
    echo(v=v, rhex=rhex);
    // proportions
    // a < 90: w1 = 2, w2 = 1, d = tan(a)/2
    // a = 90: w1 = 2, w2 = 2, d = 1/2
    // a > 90: w1 = 2, w2 = 2 - 1/tan(a), d = 1/2
    pd0 = angle < 90 ? 1/2 * tan(angle) : 1;
    pd = min(1, pd0);
    p2 = angle < 90 ? 1 + 2*(pd0-pd)/tan(angle) : 2 - 2/tan(angle);
    p1 = 2;
    a = atan2(pd, 1-p2/2);
    // fit hex to available space
    wmax = v.x - 2*ws;  // widest possible base
    whex = rhex ? min(2*rhex, wmax) : wmax;  // limit to 2*rhex
    xscale = p1 / max(p1, p2) * whex/2;
    yscale = v.y / pd;
    scale = min(xscale, yscale);
    w1 = p1 * scale;
    w2 = p2 * scale;
    d = pd * scale;
    tab([v.x, d], w1=w1, w2=w2, angle=a, rint=r, rext=r, joiner=joiner);
}
module round_tab(size=undef, d=Dthumb, r=Rext, joiner=Djoiner) {
    // approximate width of opening at the tangents
    // (quantization of $fa causes a small difference from ideal)
    axis = d/2 + r;
    span = 2*axis*cos(asin(r/axis));
    v = is_undef(size) ? area([span, d/2]) : area(size);
    intersection() {
        fillet(r, r) {
            semistadium(h=0, r=d/2);
            turnaround = [span+2*r, 2*r+EPSILON];
            translate([0, -turnaround.y/2]) square(turnaround, center=true);
        }
        translate([0, v.y/2]) square([v.x, v.y+2*joiner], center=true);
    }
}
module circle_tab(size=undef, d=Dthumb, r=Rint, joiner=Djoiner) {
    axis = d/2 + r;
    rise = d/2 - r;
    span = 2 * sqrt(axis^2 - rise^2);
    v = is_undef(size) ? area([max(d, span), d]) : area(size);
    intersection() {
        fillet(r, r) {
            translate([0, d/2]) circle(d=d);
            turnaround = [span+2*r, 2*r+EPSILON];
            translate([0, -turnaround.y/2]) square(turnaround, center=true);
        }
        translate([0, v.y/2]) square([v.x, v.y+2*joiner], center=true);
    }
}
module notch(size, w1=undef, w2=undef, angle=Atab, rint=Rint, rext=Rext, cut=Dcut) {
    // create a notch shape inside a given area
    //   size   maximum extent of notch, including base rounding
    //   w1     outer width
    //   w2     inner width (minimum)
    //   angle  rise angle
    //   rint   inner rounding
    //   rext   outer rounding
    //   cut    depth below baseline (for clean cuts)
    tab(size=size, w1=w1, w2=w2, angle=angle, rint=rext, rext=rint, joiner=cut);
}
module hex_notch(size=undef, rhex=undef, angle=Ahex, r=Rext, cut=Dcut) {
    hex_tab(size=size, rhex=rhex, angle=angle, r=r, joiner=cut);
}
module round_notch(size=undef, d=Dthumb, r=Rext, cut=Dcut) {
    round_tab(size=size, d=d, r=r, joiner=cut);
}
module circle_notch(size=undef, d=Dthumb, r=Rint, cut=Dcut) {
    circle_tab(size=size, d=d, r=r, joiner=cut);
}
module punch(d, cut=Dcut, center=false) {
    z = min(d, 0) - cut;  // allow negative cuts
    raise(z) prism(height=abs(d)+2*cut, center=center) children();
}

function wall_thickness(wall=undef, thick=false, default=Dwall) =
    let (minimum = thick ? 4*Dfpath : Dfwidth)
    max(is_undef(wall) ? pround(default) : wall, minimum);

module stacking_tabs(size, height=Htab, r=Rext, gap=Dfpath/2, slot=false) {
    v = area(size);
    h = height + Djoiner;
    d = 2*Dfpath;
    w = v.y - 2*r - 2*Dfpath;
    o = [v.x/2 - 3/2*d, w/2 - height];
    for (i=[-1,+1]) translate([o.x*i, 0]) {
        if (slot) raise(-Djoiner) hull() {
            // widen slot and slightly lengthen it
            vslot = [d+Dfpath, w+2*Dgap, h];
            prism(vslot, r=gap);
            // taper the space above the slot to ease bridging
            prism([EPSILON, vslot.y, vslot.z+2*Hflayer]);
        } else rotate([90, 0, 90]) prism(height=d, center=true)
            tab([w, height], angle=90, rint=0, rext=height);
    }
}
module scoop(size, height=undef, rint=Rint, rscoop=2*Rext, cut=Dcut) {
    v = volume(size, height);
    rmax = min(v.x/4, v.y/4, v.z);  // limit radiuses to safe values
    rn0 = min(rint, rmax);
    rn1 = min(rscoop, rmax);
    hull() {
        raise(v.z) prism(v, height=cut, r=rn0);
        for (angle=[0:$fa:90]) {
            cz = 1-cos(angle);
            cx = 1-sin(angle);
            htier = rscoop * cz;
            vtier = v - 2 * cx * area(rn1);
            rmax = min(vtier)/2 - EPSILON;
            rtier = min(rmax, cx * (rn1 - rn0) + rn0);
            raise(htier) prism(vtier, height=v.z-htier+EPSILON, r=rtier);
        }
    }
}
module stacker(size=undef, height=undef, depth=undef, r=Rext, wall=undef,
               stack=Hstack, rim=Drim, base=undef, lip=undef, gap=Dgap,
               thick=false, convex=undef) {
    v = is_undef(size) ? undef : volume(size, height);
    h = is_undef(height) ? v.z : height;
    depth = is_undef(depth) ? h - Hfloor : depth;
    wall = is_undef(wall) ? wall_thickness(wall, thick) : wall;
    rim = is_undef(rim) ? Drim : rim;
    base = stack && (is_undef(base) ? true : base);
    lip = stack && (is_undef(lip) ? true : lip);
    convex = is_undef(convex) ? !!v : convex;  // is the 2D shape convex?
    hfloor = Hfloor;
    echo(v=v, h=h, stack=stack, rim=rim, base=base, lip=lip);
    module inset(in) {
        offset(r=r-in) offset(delta=in-r)
        offset(r=-r-in) offset(delta=r)
        children();
    }
    module stack() difference() {
        if (base) {
            joiner = Dcut;  // thicker joiner to avoid clash with floors
            raise(hfloor) prism(height=h-hfloor, r=r) children();
            raise(stack) if (convex) hull() {
                prism(height=joiner) inset(rim+gap) children();
                raise(hfloor-stack) prism(height=joiner, r=r) children();
            } else minkowski() {
                prism(height=joiner) inset(rim+gap) children();
                cylinder(h=hfloor-stack, r1=0, r2=rim+gap);
            }
            prism(height=stack+joiner) inset(rim+gap) children();
        } else prism(height=h, r=r) children();
        if (lip) raise(h - stack)
            prism(height=stack+Dcut) inset(rim) children();
        if (depth) raise(h - depth)
            prism(height=depth+Dcut) inset(wall) children();
    }
    stack() if (v) square(area(v), center=true); else children();
}
module box(
    // exterior dimensions
    size=Vbox,      // exterior size in 2 or 3 dimensions (= Vbox)
    height=undef,   // wall height (overrides size.z if given)
    r=Rext,         // corner radius (= Rext)
    // wall dimensions & layout
    thick=undef,    // wall thickness selector (see below)
    wall=undef,     // exterior wall thickness (see below)
    divider=undef,  // interior wall thickness (= wall)
    grid=1,         // section layout (= [1, 1])
    // interior dimensions
    depth=undef,    // well depth (= height-Hfloor)
    scoop=false,    // scoop radius
    hole=false,     // bottom hole diameter
    // notch dimensions
    notch=false,    // top notch depth (= Dthumb)
    thumb=false,    // thumb notch width (= Dthumb)
    draw=false,     // draw notch depth (= depth)
    feet=false,     // rear feet for side draw boxes
    // stacking features
    tabs=false,     // side tabs (= size.y)
    slots=false,    // side slots (= size.y)
    stack=undef,    // stacking base height (= Hfoot)
    rim=undef,
    base=undef,
    lip=undef,
    // preview options
    color=undef,
    debug=true,
) {
    // exterior dimensions
    vbox = volume(size, height);
    height = vbox.z;
    if (debug) echo(vbox=vbox, height=height);
    // wall dimensions & layout
    thick = is_undef(thick) ? tabs || slots || thumb : thick;
    wall = is_undef(wall) ? wall_thickness(wall, thick) : wall;
    divider = wall_thickness(divider, default=wall);
    grid = grid(grid);
    if (debug) echo(thick=thick, wall=wall, divider=divider, grid=grid);
    // interior dimensions
    depth = is_undef(depth) ? vbox.z - Hfloor : depth;
    scoop = numeric_flag(scoop, 3/2*r);
    hole = numeric_flag(hole, Dthumb);
    if (debug) echo(depth=depth, scoop=scoop, hole=hole);
    vwell = volume(area(vbox) - area(2*wall), depth);
    dx = (vwell.x + divider) / grid.x;
    dy = (vwell.y + divider) / grid.y;
    vcell = vround([dx - divider, dy - divider, depth]);
    hfloor = vbox.z - vwell.z;
    if (debug) echo(vwell=vwell, vcell=vcell, hfloor=hfloor);
    // notch dimensions
    notch = numeric_flag(notch, Dthumb);
    thumb = numeric_flag(thumb, Dthumb);
    draw = numeric_flag(draw, depth);
    vcore = area(vbox) - area(4*r);  // safe cutting area
    if (debug) echo(notch=notch, thumb=thumb, draw=draw, vcore=vcore);
    // stacking features
    tabs = numeric_flag(tabs, vbox.y);
    slots = numeric_flag(slots, vbox.y);
    stack = numeric_flag(stack, Hstack);
    base = stack && (is_undef(base) ? true : base);
    lip = stack && (is_undef(lip) ? true : lip);
    if (debug) echo(tabs=tabs, slots=slots, stack=stack, base=base, lip=lip);
    // build the box
    colorize(color) difference() {
        // exterior
        union() {
            stacker(vbox, depth=0, r=r, stack=stack, base=base, lip=lip);
            if (tabs) {
                o = [0, tabs < 0 ? vbox.y/2 + tabs/2 : 0, vbox.z];
                vt = [vbox.x, abs(tabs), vbox.z];
                translate(o) stacking_tabs(vt, r=r);
            }
            if (feet) {
                o = [vcore.x/2-3/2*r, vbox.y/2+wall-r, vbox.z-3/2*r];
                for (i=[-1,+1]) scale([i, 1]) translate(o) sphere(r);
            }
        }
        // tab slots
        if (slots) {
            o = [0, slots < 0 ? vbox.y/2 + slots/2 : 0];
            vt = [vbox.x, abs(slots), vbox.z];
            translate(o) stacking_tabs(vt, r=r, slot=true);
        }
        // interior
        if (grid.x && grid.y) for (i=[1/2:grid.x]) for (j=[1/2:grid.y])
            translate([i*dx, j*dy] - area(vwell/2) - area(divider)/2) {
                if (depth) raise(hfloor) {
                    rint = r - wall;
                    if (scoop) scoop(vcell, rint=rint, rscoop=scoop, cut=Dcut);
                    else prism(vcell, height=vcell.z+Dcut, r=rint);
                }
                if (hole) raise(-Dcut) cylinder(h=hfloor+2*Dcut, d=hole);
            }
        // notch
        if (notch || draw || thumb) translate([0, -vbox.y/2, vbox.z]) {
            // long top notch (for horizontal deck boxes)
            if (notch) {
                rotate(Sdown) punch(vbox.y) hex_notch([vcore.x, notch]);
                if (lip) {
                    // fix corners on stackers
                    angle = Ahex;
                    ws = r/tan(90 - angle/2);  // shoulder width
                    wnotch = min(4*notch/tan(angle), vcore.x - 2*ws);
                    w1 = wnotch - 2*stack/tan(angle);
                    d = wnotch/2 * sin(angle) - stack;
                    raise(Djoiner-stack) rotate(Sdown) punch(vbox.y, -Drim)
                        notch([vcore.x, d], w1=w1, rint=r, rext=r);
                }
            }
            if (draw) {
                // front draw notch (for vertical deck boxes)
                angle = Adeep;
                w1 = vcore.x;
                vnotch = [vbox.x, draw];
                rotate(Sdown) punch(wall)
                    notch(vnotch, w1=w1, angle=angle, rint=r, rext=r);
                if (lip) {
                    // fix corners on stackers
                    w1 = w1 - 2*stack/tan(angle);
                    vnotch = [vnotch.x, vnotch.y-stack];
                    translate([0, Drim, Djoiner-stack])
                        rotate(Sdown) prism(height=wall-Drim+Dcut)
                        notch(vnotch, w1=w1, angle=angle, rint=r, rext=r);
                }
                if (base && depth <= draw) {
                    // soften ridge at base
                    w2 = w1 - 2*depth/tan(angle);
                    punch(-vbox.z) notch([w2, Drim+Dgap], rint=r, rext=r);
                }
            }
            if (thumb) {
                // side thumb notch (for horizontal card trays)
                punch(-vbox.z) hex_notch([vcore.x, thumb/2]);
                angle = Ahex;
                vnotch = [vcore.x, depth];
                w2 = thumb/sin(angle);
                rotate(Sdown) punch(wall) notch(vnotch, w2=w2, rext=r);
                if (lip) {
                    // fix corners on stackers
                    angle = max(angle, tab_angle(vnotch, w2, Rext));
                    vnotch = [vnotch.x, vnotch.y-stack];
                    echo(angle=angle);
                    translate([0, Drim, Djoiner-stack])
                        rotate(Sdown) prism(height=wall-Drim+Dcut)
                        notch(vnotch, w2=w2, angle=angle, rext=r);
                }
            }
        }
    }
    // children
    raise(vbox.z+EPSILON) children();
}
module box_divider(size=Vbox, height=Hcard_divider, r=Rext, wall=undef, gap=Dgap,
                   hole=false, notch=false, thick=true, color=undef) {
    vbox = area(size);
    vcore = vbox - area(2*r);
    // convert numeric flags to defaults
    hole = numeric_flag(hole, default=Dthumb);
    notch = numeric_flag(notch, default=Dthumb);
    echo(hole=hole, notch=notch);
    // dimensions
    thick = thick || notch;
    wall = wall_thickness(wall, thick);
    echo(wall=wall);
    v = area(area(vbox) - area(2*wall + 2*gap));
    echo(vbox=vbox, thick=thick, wall=wall, v=v);
    // draw box
    colorize(color) prism(height=height, r=r) difference() {
        // exterior
        square(v, center=true);
        if (hole) circle(d=hole);
        if (notch) for (i=[-1,+1]) scale([1, i])
            translate([0, -vbox.y/2]) hex_notch([vcore.x, notch/2], r=0);
    }
    // children
    raise(height+EPSILON) children();
}
module box_lid(size=Vbox, height=Hfloor, r=Rext, slots=Htab, color=undef) {
    vbox = volume(size, height);
    colorize(color) difference() {
        prism(vbox, r=r);
        stacking_tabs(vbox, height=slots, r=r, slot=true);
    }
    %stacking_tabs(vbox, height=slots, r=r);
}

module dice_rack(n=1, size=Ddice, height=undef, wall=undef, divider=undef,
                 gap=undef, color=undef) {
    vdice = volume(size);
    wall = is_undef(wall) ? pround(Dwall) : wall;
    divider = is_undef(divider) ? wall : divider;
    gap = is_undef(gap) ? ceil(wall) - wall : gap;
    cell = area(vdice) + area(2*gap);
    dx = cell.x + divider;
    w = 2*wall + dx*n - divider;
    d = 2*wall + cell.y;
    h = is_undef(height) ? ceil(vdice.z / 3) + Hfloor : height;
    box([w, d, h], grid=[n, 1], wall=wall, divider=divider, color=color);
}

module layout_tool(size, height=Hfloor, r=Rext,
                   font="Trebuchet MS:style=Bold") {
    vtool = volume(size, height, wide=true);
    label = str(size.x, "×", size.y);
    lsize = min(10, 2/3*size.y);
    echo(vtool=vtool, label=label);
    difference() {
        prism(vtool, r=r);
        raise(vtool.z - 1) prism(height=1+Dcut)
            text(label, size=lsize, font=font,
                halign="center", valign="center");
    }
}

module test_pattern() {
    v = area(Vprinter) - area(5);
    dx = v.x/2 - INCH/2;
    dy = v.y/2 - INCH/2;
    w = 4*Dfpath;
    prism(height=Hflayer) {
        for (i=[-1:1:+1]) for (j=[-1:1:+1]) translate([i*dx, j*dy]) {
            if (abs(i) == 1 && abs(j) == 1) square(INCH, center=true);
            else circle(d=INCH);
            if (i == 0) square([v.x, w], center=true);
            if (j == 0) square([w, v.x], center=true);
        }
    }
    translate([0, -dy/2]) cylinder(h=Hfloor, d=INCH);
}

module test_game_shapes() {
    // TODO: update
    %box_frame();
    module layout(i=0, j=0) {
        translate([100*i, 120*j]) children();
    }
    layout(-2) scoop_tray();
    layout(-3) grid_divider();
    layout(-4) deck_divider();
    layout(-5) tray_divider();
    layout(-6) creasing_tool();
    layout(+2) card_tray() %tray_divider();
    layout(+3) card_tray(cards=floor((Htray-Hfloor-Hlip)/Hcard));
    layout(+4) draw_tray();
    layout(+5) deck_box();
    layout(-1, +1) box([75, 60, 25], hole=true, thumb=true, stack=true);
    layout(-1, +0) box([60, 75, 25], hole=true, thumb=true, stack=true);
    layout(+0, -2) box(25, 5, grid=[2, 2], stack=true);
    layout(+0, -3/2) stacker(25, 5);
    hex3 = [[-1/2, 1], [-1/2, 0], [1/2, 0]];
    rhex = Dthumb/2;
    layout(+0, -1) hex_box(hex3, height=5, rhex=rhex, grid=rhex, hole=rhex, stack=true);
    layout(+0, +0) box([76, 50, 49], grid=2, draw=true, notch=true, stack=true);
    layout(+0, +1) box([50, 50, 49], grid=2, draw=true, notch=true, stack=true);
    layout(+0, +2) layout_tool([72, 20]);
    layout(+1, +1) box([Vtray.y, 4*Vtray.z, Vtray.x], notch=true);
    layout(+1, +0) box(volume(Vtray, wide=true), hole=true, thumb=true);
    layout(+1, -1) box([Vtray.y, Vtray.z, Vtray.x], draw=true, feet=true);
}

