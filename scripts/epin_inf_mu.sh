cd ../
EPS=0.1
SEED_DIST="uniRAND"
CAND_EDGE_PROB="orig"
SEED_MODE="IM"
K=100
for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
do
./IMAPS -func=OUTDEG -gname=Epinions -edgesize=$K -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
./IMAPS -func=PROB -gname=Epinions -edgesize=$K -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait

for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
do
    ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=0 -method="OUTDEG" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=0 -method="PROB" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait
for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
do
    ./IMAPS -func=EVAL -gname=Epinions -edgesize=0 -eps=$EPS -vanilla=1 -rand_seed=0 -method="OUTDEG" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
done
wait
for RAND_SEED in $(seq 1 5);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
    # for K in 5;
    do
        ./IMAPS -func=SANDWICH -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=RAND -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait
for RAND_SEED in $(seq 1 5);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
    # for K in 5;
    do
        ./IMAPS -func=AIS -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=SUBSIM -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 6 10);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
    # for K in 5;
    do
        ./IMAPS -func=SANDWICH -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=RAND -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 6 10);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";

    do
        ./IMAPS -func=AIS -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=SUBSIM -gname=Epinions -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 1 10);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
    do
        ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SANDWICH" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="RAND" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait

for RAND_SEED in $(seq 1 10);
do
    for SEED_DIST in "mu0.1" "mu0.3" "mu0.5" "mu0.7" "mu0.9";
    do
        ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="AIS" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
        ./IMAPS -func=EVAL -gname=Epinions -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SUBSIM" -seed_dist=$SEED_DIST -candEdges_prob=$CAND_EDGE_PROB -seed_mode=$SEED_MODE &
    done
done
wait