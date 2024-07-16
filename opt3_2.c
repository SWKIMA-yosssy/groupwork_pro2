void opt3_2(int *new_path, int n, double *new_path_length) {
  double path_length[7]; // store total length of each roots which are created
                         // by n-opt
  // i=0~2  => for 2-opt, other is for 3-opt
  int selected_nodes[3]; // store selected node numbers used in n-opt
  int opted_path[7][n];  // store total node list of each roots
  double min_path_length = *new_path_length;
  int min_path_num = -1; // number of index which gives minimum length
  // -1 means no change
  int i;

  for (i = 0; i < 3; i++) {
    selected_nodes[i] = createrand(n);
  }

  opt_2(new_path, selected_nodes, opted_path, path_length);
  opt_3(new_path, selected_nodes, opted_path, path_length);

  for (i = 0; i < 7; i++) {
    if (min_path_length > path_length[i]) {
      min_path_num = i;
      min_path_length = path_length[i];
    }
  }

  // reflect outcome
  if (min_path_num != -1) {
    *new_path_length = min_path_length;
    for (i = 0; i < n; i++) {
      new_path[i] = opted_path[min_path_num][i];
    }
  }
}
