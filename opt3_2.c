void opt_3_2_cal_length(int *pre_path, int n, int *selected_nodes,
                        double *path_length, int whether_geograph) {
  int A = selected_nodes[0], B = pre_path[(selected_nodes[0] + 1) % n];
  int C = selected_nodes[1], D = pre_path[(selected_nodes[1] + 1) % n];
  int E = selected_nodes[2], F = pre_path[(selected_nodes[2] + 1) % n];
  // initial state of node are circle: state 3 in paper
  if (whether_geograph == 0) { // data is given as euclidean data
    path_length[0] =
        euclid_distance(A, E) + euclid_distance(B, D) + euclid_distance(C, F);
    path_length[1] =
        euclid_distance(A, D) + euclid_distance(E, C) + euclid_distance(B, F);
    path_length[2] =
        euclid_distance(A, C) + euclid_distance(B, E) + euclid_distance(D, F);
    path_length[3] =
        euclid_distance(A, C) + euclid_distance(B, D) + euclid_distance(E, F);
    path_length[4] =
        euclid_distance(A, B) + euclid_distance(C, E) + euclid_distance(D, F);
    path_length[5] =
        euclid_distance(A, D) + euclid_distance(E, B) + euclid_distance(C, F);
    path_length[6] =
        euclid_distance(A, E) + euclid_distance(D, C) + euclid_distance(B, F);
  } else {
    path_length[0] = geograph_distance(A, E) + geograph_distance(B, D) +
                     geograph_distance(C, F);
    path_length[1] = geograph_distance(A, D) + geograph_distance(E, C) +
                     geograph_distance(B, F);
    path_length[2] = geograph_distance(A, C) + geograph_distance(B, E) +
                     geograph_distance(D, F);
    path_length[3] = geograph_distance(A, C) + geograph_distance(B, D) +
                     geograph_distance(E, F);
    path_length[4] = geograph_distance(A, B) + geograph_distance(C, E) +
                     geograph_distance(D, F);
    path_length[5] = geograph_distance(A, D) + geograph_distance(E, B) +
                     geograph_distance(C, F);
    path_length[6] = geograph_distance(A, E) + geograph_distance(D, C) +
                     geograph_distance(B, F);
  }
}
void opt3_2_once(int *new_path, int n, double *new_path_length, int node1,
                 int node2, int node3, int whether_geograph) {
  double path_length[7];
  int selected_nodes[3] = {node1, node2,
                           node3}; // store selected node numbers used in n-opt
                                   // int opted_path[7][n];  // store total node
                                   // list of each roots. ith global array.
  double min_path_length = *new_path_length;
  int min_path_num = -1; // number of index which gives minimum length
  // -1 means no change
  int i;

  opt_3_2_cal_length(new_path, n, selected_nodes, path_length,
                     whether_geograph);

  for (i = 0; i < 7; i++) {
    if (min_path_length > path_length[i]) {
      min_path_num = i;
      min_path_length = path_length[i];
    }
  }

  // reflect outcome
  if (min_path_num != -1) {
    *new_path_length = min_path_length;
  } else if (i == 0) {
  }
}
