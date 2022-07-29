// Game organizer library

include <cards.scad>

// TODO: normalize size and height parameters
// TODO: deck boxes from cosmic-box
// TODO: hex boxes and grids
// TODO: redesign index holes
// TODO: redesign tray feet
// TODO: parameterize uses of Vcard, Vtray, etc.

// naming conventions
// A  angle
// C  color
// D  depth, diameter, thickness
// H  height
// N  number
// P  polygon (list of points)
// Q  quality setting (draft, final)
// R  radius
// V  vector  [W, H] or [W, D, H] or [x, y, z]
// W  width

// scale adjustments:
// to counteract shrinkage, scale X & Y by 100.5% in slicer

Qdraft = 15;  // 24 segments per circle (aligns with axes)
Qfinal = 5;
$fa = Qdraft;
$fs = 0.1;

EPSILON = 0.01;
MICRON = 0.001;
PHI = (1+sqrt(5))/2;

// filament metrics
Hflayer = 0.25;
Dfwidth = 0.70;  // extrusion width
Dfoverlap = Hflayer * (1 - PI/4);  // overlap between paths
echo(Hflayer=Hflayer, Dfwidth=Dfwidth, Dfoverlap=eround(Dfoverlap, 0.001));

// organizer metrics
Dwall = 2.0;
Hfloor = Dwall;
Dgap = 0.1;
Dcut = eround(2/3*Dwall);  // cutting margin for negative spaces
echo(Dwall=Dwall, Hfloor=Hfloor, Dgap=Dgap, Dcut=Dcut);
Rext = 3.0;  // external corner radius
Rint = Rext - Dwall;  // internal corner radius
echo(Rext=Rext, Rint=Rint);
Avee = 60;  // angle for index v-notches
Dthumb = 25.0;  // index hole diameter
Dstrut = 12.0;  // width of struts and corner braces (TODO: rework or remove)
echo(Avee=Avee, Dthumb=Dthumb, Dstrut=Dstrut);

// game box interior
Vgame = [288, 288, 69];  // typical FFG box interior
Hwrap = 55;  // cover art wrap ends here, approximately
echo(Vgame=Vgame, Hwrap=Hwrap);
Hmanual = 1.0;
Hceiling = Vgame.z - eceil(Hmanual, 0.5);
echo(Hmanual=Hmanual, Hceiling=Hceiling);

// card metrics
Vcard = [66, 91];  // Gamegenic gray sleeves
Hcard = 0.525;  // generic card with Gamegenic prime sleeves
echo(Vcard=Vcard, Hcard=Hcard);

// container metrics
echo(card_pile_height(60, Hcard_universal)+Hfloor+2*Hcard_divider);  // TODO
Htray = 13;
Vtray = [72, 100, Htray];
Hdeck = 65.5;  // TODO: remove
Vlongtray = [72, 144, Htray];  // TODO: remove
// nesting feet (TODO: reshape, maybe to area(Vtray)/10)
Rfoot = Rint - Dgap;  // concentric with Rint & Rext with nesting gap
Dfoot = 10.0 - 2*Dgap;
Hfoot = 1.0;
echo(Rfoot=Rfoot, Dfoot=Dfoot, Hfoot=Hfoot);

// TODO: eliminate these
Dlong = Vgame.y / 2;
Vlong = deck_box_volume(Dlong);
Dshort = Vgame.x - 4*Vlong.x;
echo(cards=floor(Dlong / Hcard));


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

// utility functions
function sum(v) = v ? [for(p=v) 1]*v : 0;
function swapxy(v) = [v.y, v.x, if (2<len(v)) for (i=[2:len(v)-1]) v[i]];
function unit_axis(n) = [for (i=[0:2]) i==n ? 1 : 0];

// TODO: review these
function deck_box_volume(d) = [Vcard.y + 2*Rext, d, Hdeck];
function tray_volume(height=1) =
    [Vtray.x, Vtray.y, height*Vtray.z];
function deck_height(n=0, deck=1, gap=Dgap/2) =
    deck*(Hdeck+gap) + n*(Htray+gap);
function card_pile_height(n, height=Hcard) = n*height;

// transformations
module colorize(c=undef, alpha=undef) {
    // skip the color() call if both parameters are undef
    if (is_undef(c) && is_undef(alpha)) children();
    else color(c, alpha) children();
}
module flatten(size, space) {
    // shear and flatten with fixed sides (like flattening a cardboard box)
    // TODO: handle low walls
    v = volume(size);
    dx = min(space - v.x, v.z);
    a = acos(dx/v.z);
    mlean = [
        [1, 0, cos(a), 0],
        [0, 1, 0, 0],
        [0, 0, sin(a), 0],
    ];
    multmatrix(m=mlean) children();
}
module lean(size, space, a, enforce=false) {
    // shear and rotate with fixed volume (like leaning cards against a box)
    // TODO: handle low walls
    // TODO: solve for a
    v = volume(size);
    A = 90 -a;
    x = v.x/sin(A);
    c = x + v.z*cos(A);
    mshear = [
        [1, 0, 0, v.x/2],
        [0, 1, 0, 0],
        [tan(a), 0, 1, v.x*tan(a)/2],
    ];
    mrotate = [
        [cos(a), 0, sin(a), -v.x/2],
        [0, 1, 0, 0],
        [-sin(a), 0, cos(a), 0],
    ];
    multmatrix(m=mrotate)
    multmatrix(m=mshear)
    children();
    // diagnostics: show how close the lean is to ideal
    echo(a=a, footprint=mround(x), lean=mround(c), gap=mround(space-c));
    if (enforce) assert(mround(space-c) == 0);
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
module stadium(h, r=undef, d=undef) {
    // creates a stadium with rectangle width w and radius r,
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
// TODO: capsule and semicapsule modules with cylinder height h and radius r

module prism(size=undef, height=undef, r=undef, rint=undef, rext=undef,
             scale=1, center=false) {
    v = is_undef(size) ? undef : volume(size, height);
    h = is_undef(height) ? size.z : height;
    ri = is_undef(rint) ? is_undef(r) ? 0 : r : rint;  // inside turns
    re = is_undef(rext) ? is_undef(r) ? 0 : r : rext;  // outside turns
    linear_extrude(height=h, scale=scale, center=center) fillet(ri, re) {
        if (is_undef(v)) children();
        else square(area(v), center=true);
    }
}
module box_frame(size=Vgame, wall=Dwall, wrap=Hwrap, gap=Dgap) {
    // create the outline of a box with given interior and thickness
    vint = volume(size);
    vext = vint + [2*wall, 2*wall, wall];
    dwall = wall - gap/2;  // shrink the wall to leave a small gap
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

module wall_vee_cut(size, a=Avee, cut=Dcut) {
    // TODO: round bottom of vee
    span = size.x;
    y0 = -2*Rext;
    y1 = size.z;
    rise = y1;
    run = a == 90 ? 0 : rise/tan(a);
    x0 = span/2;
    x1 = x0 + run;
    a1 = (180-a)/2;
    x2 = x1 + Rext/tan(a1);
    x3 = x2 + Rext + EPSILON;  // needs +EPSILON for 90-degree angles
    poly = [[x3, y0], [x3, y1], [x1, y1], [x0, 0], [x0, y0]];
    rotate([90, 0, 0]) linear_extrude(size.y+2*cut, center=true)
    difference() {
        translate([0, y1/2+cut/2]) square([2*x2, y1+cut], center=true);
        for (s=[-1,+1]) scale([s, 1]) hull() {
            fillet(Rint, Rext) polygon(poly);
            translate([x0, y0]) square([x3-x0, -y0]);
        }
    }
}

module deck_box(d=Dlong, seed=undef, tiers=2, flip=false, color=undef) {
    // TODO: add optional draw cut and tilt feet
    vbox = deck_box_volume(d);
    shell = area(vbox);
    well = shell - 2 * area(Dwall);
    hole = shell - 2 * area(Dstrut);
    // notch dimensions:
    dtop = vbox.x - 2*Dstrut;  // corner supports
    hvee = vbox.z - dtop/2 * sin(Avee);
    vend = [dtop/2, vbox.y, vbox.z-hvee];
    colorize(color) difference() {
        // outer shell
        prism(vbox, r=Rext);
        // card well
        raise(Hfloor) prism(well, height=vbox.z, r=Rint);
        // base round
        raise(-Dgap) prism(hole, height=vbox.z, r=Dthumb/2);
        // side cuts
        raise(hvee) wall_vee_cut(vend);  // end vee
    }
}
module creasing_tool(cards=10) {
    // block for creasing index wrappers
    margin = (Vcard.x - Dthumb) / 2;
    linear_extrude(height=card_pile_height(cards)) difference() {
        square(Vcard, center=true);
        stadium_fill([Dthumb, Vcard.y - 2*margin]);
    }
}

module tray_feet_cut(tray=Vtray, cut=Dcut) {
    // d = distance from center of foot to edge of tray
    // (half width of foot + external radius - concentric foot radius)
    d = Dfoot/2 + Rext - Rfoot;
    for (i=[-1,+1]) for (j=[-1,+1])
        translate([i * (tray.x/2-d), j * (tray.y/2-d)])
            tray_foot(cut=cut);
}
module tray_foot(cut=0) {
    hleg = Hfloor - Hfoot;
    dleg = Dfoot - 2*Rfoot;
    if (cut) {
        raise(-cut) prism(dleg, height=hleg+cut+Hflayer);
        %raise(-Hfoot) tray_foot();
    } else {
        prism(Dfoot, height=Hfoot, r=Rfoot);
        prism(dleg, height=Hfoot+hleg);
    }
}
module card_well(h=1, slope=false, cut=Dcut) {
    // TODO: remove slope option
    vtray = tray_volume(height=h);
    shell = area(vtray);
    well = shell - 2 * area(Dwall);
    raise(Hfloor) {
        hwell = vtray.z - Hfloor;
        h1 = hwell - Hfoot;
        v1 = well - 2 * area(Rint);
        // if set, the well slopes inward to center the cards
        v0 = slope ? Vcard : v1;
        hull() {
            raise(hwell) prism(well, height=cut, r=Rint);
            minkowski() {
                sphere(Rint);
                hull() {
                    raise(hwell-Rint) prism(v1, height=Rint);
                    raise(Rint) prism(v0, height=h1-2*Rint);
                }
            }
        }
        // thumb vee
        dy = vtray.z - Hfloor;
        dx = vtray.x/2 - Dthumb/2 - Dstrut/2;
        a = max(Avee, atan(dy/dx));
        translate([0, Dwall-vtray.y]/2)
            wall_vee_cut([Dthumb, Dwall, vtray.z-Hfloor], a=a, cut=cut);
    }
    raise(-cut) linear_extrude(vtray.z+2*cut-EPSILON) {
        // thumb round
        xthumb = 2/3 * Dthumb;  // depth of thumb round
        translate([0, -cut-vtray.y/2])
            semistadium(xthumb-Dthumb/2+cut, d=Dthumb);
        // bottom index hole
        translate([0, xthumb/3])
        stadium(xthumb, d=Dthumb);
    }
}
module card_tray(h=1, cards=0, color=undef) {
    hfloor = max(Hfloor, 2.0*h/2);
    vtray = tray_volume(height=h);
    shell = [vtray.x, vtray.y];
    well = shell - 2 * area(Dwall);
    origin = [well.x/2 + Dwall - shell.x/2, 0];
    dx = well.x + Dwall;
    colorize(color) difference() {
        prism(vtray, r=Rext);
        card_well(h=h);
        tray_feet_cut();
    }

    %raise() // card stack
        if (cards) card_pile(cards) children();
        else children();
}
module draw_tray(h=2, cards=0, color=undef) {
    vtray = tray_volume(height=h+1);
    shell = [vtray.x, vtray.y];
    well = shell - 2 * area(Dwall);
    origin = [well.x/2 + Dwall - shell.x/2, 0];
    dx = well.x + Dwall;
    srise = Htray - Hfloor;
    slope = atan(srise / well.y);
    mslope = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, -sin(slope), 1, 0],
    ];
    colorize(color) difference() {
        prism(height=vtray.z, r=Rext) difference() {
            square(shell, center=true);
            xthumb = 2/3 * Dthumb;
            // thumb round
            translate([0, -Dcut - Vtray.y/2])
                semistadium(xthumb-Dthumb/2+Dcut, d=Dthumb);
            // bottom index hole
            translate([0, xthumb/3])
            stadium(xthumb, d=Dthumb);
        }
        tray_feet_cut();
        // sloped floor
        raise(Hfloor + srise/2) multmatrix(m=mslope)
            prism(well, height=vtray.z+Dcut, r=Rint);
        // open front
        translate([0, -vtray.y/2, Htray]) rotate(90)
            wall_vee_cut([2*Rext, vtray.x + Dcut, vtray.z-Htray]);
    }
}
module card_pile(n=10, up=false, color=undef) {
    // TODO: generalize this & don't repeat so much code
    hcards = card_pile_height(n);
    vcard = area(Vcard);
    spin = up ? [0, 90, 0] : 0;
    lift = up ? vcard.x/2 : 0;
    raise(lift) rotate(spin) colorize(color) prism(vcard, height=hcards);
    translate(spin ? [hcards, 0, 0] : [0, 0, hcards]) children();
}
module deck_divider(v=Vcard_divider, h=Hcard_divider, up=false, color=undef) {
    // vertical divider for wide deck boxes
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
        prism(height=h, r=Rext) polygon(poly);
    }
    translate(spin ? [h, 0, 0] : [0, 0, h]) children();
}
module tray_divider(v=Vcard_divider, h=Hcard_divider, color=undef) {
    // horizontal divider for card trays
    colorize(color) prism(height=h, r=Rext) difference() {
        square(v, center=true);
        xthumb = 2/3 * Dthumb;  // depth of thumb round
        // match the index holes in the underlying card tray
        // thumb round
        translate([0, -Dcut - Vtray.y/2])
            semistadium(xthumb-Dthumb/2+Dcut, d=Dthumb);
        // bottom index hole
        translate([0, xthumb/3])
        stadium(xthumb, d=Dthumb);
    }
    translate([0, 0, h]) children();
}
module scoop_well(v, h, r0, r1, lip=Hfoot+Dgap, cut=Dcut) {
    hmax = h - lip;  // leave room for nesting feet
    rmax = min(v) / 4;  // limit radiuses to safe values
    rn0 = min(r0, rmax);
    rn1 = min(r1, rmax);
    hull() {
        raise(h) prism(v, height=cut, r=rn0);
        for (a=[0:$fa:90]) {
            cz = 1-cos(a);
            cx = 1-sin(a);
            htier = hmax * cz;
            vtier = v - 2 * cx * area(rn1);
            rmax = min(vtier)/2 - EPSILON;
            rtier = min(rmax, cx * (rn1 - rn0) + rn0);
            raise(htier) prism(vtier, height=h-htier+EPSILON, r=rtier);
        }
    }
}
module token_tray(scoop=2*Rext, color=undef) {
    // TODO: configurable partitions
    vtray = Vtray;
    shell = [vtray.x, vtray.y];
    origin = [Dwall-shell.x/2, Dwall-shell.y/2];
    wella = [vtray.x-2*Dwall, Vlongtray.y - vtray.y - Dwall];
    wellb = [(vtray.x-3*Dwall)/2, vtray.y - wella.y - 3*Dwall];
    colorize(color) difference() {
        prism(vtray, r=Rext);
        raise(Hfloor) {
            walls = 2 * area(Dwall);
            dva = (vtray - wella - walls) / 2;
            dvb = (wellb - vtray + walls) / 2;
            translate(dva)
                scoop_well(wella, vtray.z-Hfloor, r0=Rint, r1=scoop);
            for (i=[-1,+1]) translate([i*dvb.x, dvb.y])
                scoop_well(wellb, vtray.z-Hfloor, r0=Rint, r1=scoop);
        }
        tray_feet_cut();
    }
    %raise() rotate(90) children();  // card stack
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

module test_game_shapes() {
    %box_frame();
    module grid(i=0, j=0) {
        translate([100*i, 100*j]) children();
    }
    grid(+0) deck_box(Dlong);
    grid(-1) token_tray();
    grid(-2) deck_divider();
    grid(-3) tray_divider();
    grid(-4) creasing_tool();
    grid(+1) card_tray(h=1, cards=15);
    grid(+2) card_tray(h=2, cards=30);
    grid(+3) draw_tray(h=2, cards=30);
    grid(+0, -1.25) layout_tool([72, 20]);
}

