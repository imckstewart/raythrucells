#!/usr/bin/perl -w

use strict;

my %module_names;
my $status;

push @{$module_names{'rt_utils'}}, 'gramSchmidt';
push @{$module_names{'rt_utils'}}, 'calcDotProduct';
push @{$module_names{'rt_utils'}}, 'calcCellCentres';
push @{$module_names{'rt_utils'}}, '_getNextEdgeSet';
push @{$module_names{'rt_utils'}}, '_gridPointsAreNeighbours';
push @{$module_names{'rt_utils'}}, '_edgesFormACell';
push @{$module_names{'rt_utils'}}, '_cellVerticesMatch';
push @{$module_names{'rt_utils'}}, '_addRawCell';
push @{$module_names{'rt_utils'}}, 'getCellsFromGrid';
push @{$module_names{'rt_utils'}}, 'getEdges';
push @{$module_names{'rt_utils'}}, 'calcBaryCoords';

push @{$module_names{'raythrucells'}}, '_extractFace';
push @{$module_names{'raythrucells'}}, '_getNewEntryFaceI';
push @{$module_names{'raythrucells'}}, '_calcFaceInNMinus1';
push @{$module_names{'raythrucells'}}, '_intersectLineWithFace';
push @{$module_names{'raythrucells'}}, '_followGoodChain';
push @{$module_names{'raythrucells'}}, '_buildRayCellChain';
push @{$module_names{'raythrucells'}}, 'followRayThroughCells';

push @{$module_names{'meshtocube'}}, '_interpolateAtFace';
push @{$module_names{'meshtocube'}}, '_getFaceInterpsAlongRay';
push @{$module_names{'meshtocube'}}, '_interpOnGridAlongRay';
push @{$module_names{'meshtocube'}}, '_generateVoxelIndex';
push @{$module_names{'meshtocube'}}, '_generateNextPixelCombo';
#push @{$module_names{'meshtocube'}}, 'cellsToHyperCube';

push @{$module_names{'second_order'}}, '_evaluate2ndOrderShapeFns';
push @{$module_names{'second_order'}}, '_interpolate2ndOrderCell';
push @{$module_names{'second_order'}}, '_getParabolicShapeFns';
push @{$module_names{'second_order'}}, 'interpolateParabolic';
push @{$module_names{'second_order'}}, '_faceBaryToCellBary';
push @{$module_names{'second_order'}}, '_fillBaryBuffValues';
push @{$module_names{'second_order'}}, '_getInterpsAlongRay';
#push @{$module_names{'second_order'}}, '_setRasterFlags';
push @{$module_names{'second_order'}}, 'interpOnGridAlongRay2ndOrder';

OUTER: for my $module_name (keys %module_names){
  for my $func_name (@{$module_names{$module_name}}){
#    my $command = "../mytest $module_name $func_name";
    my $command = "./mytest $module_name $func_name";
    print "$command\n";
    $status = system($command);
    if($status){
      print "  Test failed, exit status=$status\n";
last OUTER;
    }else{
      print "  Test passed.\n\n";
    }
  }
}
if(!$status){
  print "*** Passed all tests. ***\n"
}

