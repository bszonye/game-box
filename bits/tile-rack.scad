echo("\n\n====== TILE RACK ======\n\n");

include <game-box/game.scad>

Qprint = Qfinal;  // or Qdraft

// Bloodstones: 4-tile rack
tile_rack(n=4, size=[40, 20, 6.7], $fa=Qprint);
