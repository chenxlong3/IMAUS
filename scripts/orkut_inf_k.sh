cd ../
EPS=0.1
SEED_DIST="uniRAND"
CAND_EDGE_PROB="orig"
SEED_MODE="IM"

for K in 5 10 20 40 60 80 100;
do
./IMAPS -func=OUTDEG -gname=orkut -edgesize=$K -eps=$EPS -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
./IMAPS -func=PROB -gname=orkut -edgesize=$K -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait

for K in 5 10 20 40 60 80 100;
do
    ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=0 -method="OUTDEG" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=0 -method="PROB" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait
for K in 0;
do
    ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=0 -method="OUTDEG" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait
for RAND_SEED in $(seq 1 5);
do
    for K in 5 10 20 40 60 80 100;
    # for K in 5;
    do
        ./IMAPS -func=SANDWICH -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=RAND -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait
for RAND_SEED in $(seq 1 5);
do
    for K in 5 10 20 40 60 80 100;
    do
        ./IMAPS -func=AIS -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=SUBSIM -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 6 10);
do
    for K in 5 10 20 40 60 80 100;
    # for K in 5;
    do
        ./IMAPS -func=SANDWICH -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=RAND -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 6 10);
do
    for K in 5 10 20 40 60 80 100;
    do
        ./IMAPS -func=AIS -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=SUBSIM -gname=orkut -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 1 10);
do
    for K in 5 10 20 40 60 80 100;
    do
        ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SANDWICH" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="RAND" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 1 10);
do
    for K in 5 10 20 40 60 80 100;
    do
        ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="AIS" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=EVAL -gname=orkut -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SUBSIM" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait