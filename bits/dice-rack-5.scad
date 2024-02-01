echo("\n\n====== DICE RACK (5) ======\n\n");

include <game-box/game.scad>

Qprint = Qfinal;  // or Qdraft

// Yahtzee rack, sized for Altoids tin
dice_rack(n=5, $fa=Qprint);
*dice_rack(n=5, div=1.25, $fa=Qprint);
