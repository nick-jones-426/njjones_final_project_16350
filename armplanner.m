function[armplan, armplanlength, nodesgenerated, plantime] = armplanner(envmap, armstart, armgoal, planner_id)
%call the planner in C
[armplan, armplanlength, nodesgenerated, plantime] = planner(envmap, armstart, armgoal, planner_id);

