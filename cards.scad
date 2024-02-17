// Common card and sleeve sizes

INCH = 25.4;

// Standard (poker size)
Vcard_poker = [2.5*INCH, 3.5*INCH];
Hcard_poker = 0.300;  // Bicycle deck measures ca. 16.0 mm / 54 cards
Hcard_radlands = 0.375;  // measured ca. 375 microns
// American Standard (bridge/whist size)
Vcard_bridge = [2.25*INCH, 3.5*INCH];
Hcard_bridge = 0.300;  // Bicycle deck measures ca. 16.0 mm / 54 cards
Hcard_cosmic = 0.325;  // measured ca. 325 microns
// Tarot
Vcard_tarot = [2.75*INCH, 4.75*INCH];
// French Tarot
Vcard_french = [61, 112];  // Fournier deck measures 60 x 110 mm
Hcard_french = 0.300;  // Fournier deck measures ca. 24.0 mm / 78 cards
// Standard European
Vcard_euro = [59, 92];  // e.g. Ravensburger playing cards
Hcard_dominion = 0.320;  // measured ca. 320 microns
Hcard_targi = 0.340;  // measured ca. 540 microns with Gamegenic sleeves
// Cosmic Encounter alien sheets
Vcard_alien = [115, 171];
Hcard_alien = 0.365;  // measured 355-370 microns
// other
Hcard_index = 0.25;  // cardstock for dividers & wrappers (approximate)

// sleeve metrics
// Gamegenic
Vsleeve_sand = [81, 122];  // Dixit
Vsleeve_orange = [73, 122];  // Tarot
Vsleeve_magenta = [72, 112];  // Scythe
Vsleeve_brown = [67, 103];  // 7 Wonders
Vsleeve_lime = [82, 82];  // Big Square
Vsleeve_blue = [73, 73];  // Square
Vsleeve_dark_blue = [53, 53];  // Mini Square
Vsleeve_gray = [66, 91];  // Standard Card
Vsleeve_purple = [62, 94];  // Standard European
Vsleeve_ruby = [46, 71];  // Mini European
Vsleeve_green = [59, 91];  // Standard American
Vsleeve_yellow = [44, 67];  // Mini American
Vsleeve_catan = [56, 82];  // Catan (English)
// Sleeve Kings
Vsleeve_card_game = [66.5, 91];
Vsleeve_euro = [61.5, 94];  // Standard European
Vsleeve_mini_euro = [47, 70];  // Mini European
Vsleeve_super_large = [104, 129];
// thickness
Hsleeve_none = 0;
Hsleeve_penny = 0.08;  // 40 micron sleeves (Mayday standard, UG classic)
Hsleeve_inner = 0.10;  // 50 micron sleeves (inner sleeves, UG premium soft)
Hsleeve_kings = 0.13;  // 60 micron sleeves (Sleeve Kings standard)
Hsleeve_premium = 0.18;  // 90 micron sleeves (Mayday premium)
Hsleeve_prime = 0.20;  // 100 micron sleeves (Gamegenic prime)
Hsleeve_supreme = 0.23;  // 115 micron sleeves (Ultimate Guard supreme)
Hsleeve_dragon = 0.24;  // 120 micron sleeves (Dragon Shield)
Hsleeve_double = 0.30;  // 100 + 50 micron double sleeve

// bounding box for common sleeved cards
Vcard_universal = [66, 94];  // fits sleeved standard CG & euro cards

// deck box metrics (approximate inner dimensions)
Vomnihive = [399, 106, 82];  // x2
Homnihive_rail = 57;  // height of center & side rails
Vomnihive_tray = [72, 96, 36];
Homnihive_tray_notch = 23;  // height of thumb notch
Varkhive = [303, 76.5, 102];
Varkhive_rail = 43;  // height of side rails

// card metrics: override these!
Hcard_unsleeved = 0.325;  // typical unsleeved thickness
Hcard_sleeve = Hsleeve_prime;
Vcard = Vcard_universal;
Hcard = Hcard_unsleeved + Hcard_sleeve;
echo(Vcard=Vcard, Hcard=Hcard);

// card dividers: override these!
Vcard_divider = [67, 93];
Hcard_divider = 1.0;
