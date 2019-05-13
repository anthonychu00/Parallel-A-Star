# Parallel-A-

The following is a parallelized implementation of the A* algorithm, which is used in many modern day applications such as pathfinding and protein design. As opposed to a single priority queue in the normal A* algorithm, we use multiple threads which work on different priority queues. Thus you can expand explored states at a much faster rate.
