cd ../
EPS=0.1
RAND_SEED=0
K=100
SEED_MODE="IM"
# ./IMAUS -func=format -gname=Epinions -seedsize=50 -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED

# ./IMAUS -func=SANDWICH -gname=facebook -seedsize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED
# ./IMAUS -func=SANDWICH -gname=DBLP -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED
# ./IMAUS -func=SANDWICH -gname=DBLP -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED
# ./IMAUS -func=GREEDY -gname=facebook -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -num_samples=1000 -seed_mode="IM"
# ./IMAUS -func=OUTDEG -gname=facebook -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_mode="RAND"
# ./IMAUS -func=EVAL -gname=facebook -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="GREEDY" -seed_mode="IM"

# ./IMAUS -func=OUTDEG -gname=EC -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -num_samples=1000 -seed_mode="RAND"
# ./IMAUS -func=EVAL -gname=EC -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="OUTDEG" -seed_mode="RAND"
# ./IMAUS -func=GREEDY -gname=EC -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -num_samples=10000 -seed_mode=$SEED_MODE
# ./IMAUS -func=EVAL -gname=EC -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="GREEDY" -seed_mode=$SEED_MODE
./IMAUS -func=SANDWICH -gname=facebook -edgesize=100 -eps=$EPS -rand_seed=0 -seed_dist=uniRAND -seed_mode=$SEED_MODE
# ./IMAUS -func=OUTDEG -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -seed_dist=uniRAND
# ./IMAUS -func=PROB -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -seed_dist=uniRAND
# ./IMAUS -func=AIS -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -seed_dist=uniRAND
# ./IMAUS -func=RAND -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -seed_dist=uniRAND
# ./IMAUS -func=SUBSIM -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="SANDWICH" -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="OUTDEG" -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="PROB" -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="AIS" -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="RAND" -seed_dist=uniRAND
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=100 -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE -method="SUBSIM" -seed_dist=uniRAND

# ./IMAUS -func=SUBSIM -gname=DBLP -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE &
# ./IMAUS -func=RAND -gname=DBLP -edgesize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE &

# ./IMAUS -func=RAND -gname=EC -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -num_samples=50000 -seed_mode="RAND"
# ./IMAUS -func=EVAL -gname=EC -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="RAND" -seed_mode="RAND"

# ./IMAUS -func=SANDWICH -gname=EC -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -num_samples=50000 -seed_mode="IM"
# ./IMAUS -func=EVAL -gname=EC -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SANDWICH" -seed_mode="IM"

# ./IMAUS -func=EVAL -gname=facebook -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="GREEDY" -seed_mode="RAND"
# ./IMAUS -func=EVAL -gname=DBLP -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="SUBSIM"
# ./IMAUS -func=EVAL -gname=facebook -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="GREEDY" -seed_mode="RAND"

# ./IMAUS -func=RAND -gname=facebook -seedsize=$K -eps=$EPS -rand_seed=$RAND_SEED -seed_mode=$SEED_MODE
# ./IMAUS -func=EVAL -gname=facebook -edgesize=$K -eps=$EPS -vanilla=1 -rand_seed=$RAND_SEED -method="RAND" -seed_mode=$SEED_MODE
