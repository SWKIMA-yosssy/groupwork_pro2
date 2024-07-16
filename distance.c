double distance(
    int *path, int n,
    int whether_geograph) { // int *path is array of root, n is number of node
  // int whether_geograph = 0 means data is given as euclidean, =1 means as
  // geographical
  double total_distance = 0;
  int i;

  if (whether_geograph ==
      0) { // calculate length of whole root of euclidean data
    for (i = 0; i < n - 1; i++) {
      total_distance += euclid_distance(path[i], path[i + 1]);
    }
    total_distance += euclid_distance(path[n - 1], path[0]);
  } else {
    for (i = 0; i < n - 1; i++) {
      total_distance += geograph_distance(path[i], path[i + 1]);
    }
    total_distance += geograph_distance(path[n - 1], path[0]);
  }
  return total_distance;
}
