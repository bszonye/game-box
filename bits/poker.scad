echo("\n\n====== POKER ACCESSORIES ======\n\n");

include <game-box/game.scad>

Qprint = Qfinal;  // or Qdraft

Hfloor = 1.25;  // adds up to 18mm with 40mm chips

chip_tray(n=25, rows=2, $fa=Qprint);
*chip_tray(n=5, rows=1, $fa=Qprint);
*chip_tray($fa=Qprint);
