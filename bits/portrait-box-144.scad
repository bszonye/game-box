echo("\n\n====== PORTRAIT DECK BOX ======\n\n");

include <game-box/game.scad>

Qprint = Qfinal;  // or Qdraft

Vcard_divider = [93, 67];
Hcard_divider = 1.0;
Hfoot = 0;

deck_box(size=[93, 66], width=144, $fa=Qprint);
*deck_divider($fa=Qprint);
