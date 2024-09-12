echo("\n\n====== TEST PATTERN ======\n\n");

include <game-box/game.scad>

Qprint = Qfinal;  // or Qdraft

test_pattern($fa=Qprint);
