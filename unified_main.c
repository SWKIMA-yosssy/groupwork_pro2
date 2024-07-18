#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define maxN 500
#define inf 1000000
#define earthr 6378.388  // 地球の半径
double dist[maxN][maxN]; // 重み行列
                         // visited[maxN]

int createrand(int n) { // n is number of nodes
  return rand() % n;
}
/*ユークリッド*/
struct city {
  double x;
  double y;
};

double euclid_distance(struct city *cities, int a, int b) {
  double dx = cities[a].x - cities[b].x;
  double dy = cities[a].y - cities[b].y;
  return sqrt(dx * dx + dy * dy);
}

/*地理的距離*/

// 度をラジアンに変換する関数
double deg2rad(double deg) { return deg * (M_PI / 180.0); }

/*地理的距離を求める関数*/
double geograph_distance(struct city *cities, int a, int b) {
  double dlat = deg2rad(cities[a].x - cities[b].x);
  double dlon = deg2rad(cities[a].y - cities[b].y);
  double A = sin(dlat / 2) * sin(dlat / 2) + cos(deg2rad(cities[b].x)) *
                                                 cos(deg2rad(cities[a].x)) *
                                                 sin(dlon / 2) * sin(dlon / 2);
  double C = 2 * atan2(sqrt(A), sqrt(1 - A));
  return earthr * C;
}
double distance(
    struct city *cities, int *path, int n,
    int whether_geograph) { // int *path is array of root, n is number of node
  // int whether_geograph = 0 means data is given as euclidean, =1 means as
  // geographical
  double total_distance = 0;
  int i;

  if (whether_geograph ==
      0) { // calculate length of whole root of euclidean data
    for (i = 0; i < n - 1; i++) {
      total_distance += euclid_distance(cities, path[i], path[i + 1]);
    }
    total_distance += euclid_distance(cities, path[n - 1], path[0]);
  } else {
    for (i = 0; i < n - 1; i++) {
      total_distance += geograph_distance(cities, path[i], path[i + 1]);
    }
    total_distance += geograph_distance(cities, path[n - 1], path[0]);
  }
  return total_distance;
}
/*mainで必要なもの*/

/*distとvisitedの初期化*/
// for(int i = 0;i < N ;i++){for(int j;j<N;j++){dist[i][j] == -1}}
// for (int i = 0; i < N; i++){ visited[i] = -1}

/*呼び出し*/
// while (path_size){int farthest = find_farthest_node(dist,N<visited);}
// insert_node(new_path,&path_size,dist,farthest);
// visited[farthest]=1;}

// 最も遠い点を見つける関数
int find_farthest_node(int N, int visited[maxN]) {
  int farthest = -1;
  double max_dist = -1;
  for (int i = 0; i < N; i++) {
    if (visited[i] == -1) {
      for (int j = 0; j < N; j++) {
        if (visited[j] == -1 && dist[i][j] > max_dist) {
          max_dist = dist[i][j];
          farthest = i;
        }
      }
    }
  }
  return farthest;
}

// 点をツアーに挿入する関数
void insert_nodes(int new_path[maxN], int *path_size, int node) {
  int best_position = -1;
  double best_cost;

  for (int i = 0; i < *path_size; i++) {
    int j = (i + 1) % *path_size;
    double cost = dist[new_path[i]][node] + dist[node][new_path[j]] -
                  dist[new_path[i]][new_path[j]];
    if (cost < best_cost) {
      best_cost = cost;
      best_position = j;
    }
  }

  for (int i = *path_size; i > best_position; i--) {
    new_path[i] = new_path[i - 1];
  }
  new_path[best_position] = node;
  (*path_size)++;
}
void opt_3_2_cal_length(struct city *cities, int *pre_path, int n,
                        int *selected_nodes, double *path_length,
                        int whether_geograph) {
  int A = pre_path[selected_nodes[0]],
      B = pre_path[(selected_nodes[0] + 1) % n];
  int C = pre_path[selected_nodes[1]],
      D = pre_path[(selected_nodes[1] + 1) % n];
  int E = pre_path[selected_nodes[2]],
      F = pre_path[(selected_nodes[2] + 1) % n];
  // initial state of node are circle: state 3 in paper
  if (whether_geograph == 0) { // data is given as euclidean data
    path_length[0] = euclid_distance(cities, A, E) +
                     euclid_distance(cities, B, D) +
                     euclid_distance(cities, C, F);
    path_length[1] = euclid_distance(cities, A, D) +
                     euclid_distance(cities, E, C) +
                     euclid_distance(cities, B, F);
    path_length[2] = euclid_distance(cities, A, C) +
                     euclid_distance(cities, B, E) +
                     euclid_distance(cities, D, F);
    path_length[3] = euclid_distance(cities, A, C) +
                     euclid_distance(cities, B, D) +
                     euclid_distance(cities, E, F);
    path_length[4] = euclid_distance(cities, A, B) +
                     euclid_distance(cities, C, E) +
                     euclid_distance(cities, D, F);
    path_length[5] = euclid_distance(cities, A, D) +
                     euclid_distance(cities, E, B) +
                     euclid_distance(cities, C, F);
    path_length[6] = euclid_distance(cities, A, E) +
                     euclid_distance(cities, D, C) +
                     euclid_distance(cities, B, F);
  } else {
    path_length[0] = geograph_distance(cities, A, E) +
                     geograph_distance(cities, B, D) +
                     geograph_distance(cities, C, F);
    path_length[1] = geograph_distance(cities, A, D) +
                     geograph_distance(cities, E, C) +
                     geograph_distance(cities, B, F);
    path_length[2] = geograph_distance(cities, A, C) +
                     geograph_distance(cities, B, E) +
                     geograph_distance(cities, D, F);
    path_length[3] = geograph_distance(cities, A, C) +
                     geograph_distance(cities, B, D) +
                     geograph_distance(cities, E, F);
    path_length[4] = geograph_distance(cities, A, B) +
                     geograph_distance(cities, C, E) +
                     geograph_distance(cities, D, F);
    path_length[5] = geograph_distance(cities, A, D) +
                     geograph_distance(cities, E, B) +
                     geograph_distance(cities, C, F);
    path_length[6] = geograph_distance(cities, A, E) +
                     geograph_distance(cities, D, C) +
                     geograph_distance(cities, B, F);
  }
}

void reverse(int *new_path, int n, int start, int end) {
  int temp;
  while (start < end) {
    temp = new_path[start];
    new_path[start] = new_path[end];
    new_path[end] = temp;
    if (start == n - 1) {
      start = 0;
    } else {
      start++;
    }
    if (end == 0) {
      end = n - 1;
    } else {
      end--;
    }
  }
}
int opt3_2_once(struct city *cities, int *new_path, int n,
                double *new_path_length, int node1, int node2, int node3,
                int whether_geograph) {
  double path_length[7];
  int selected_nodes[3] = {node1, node2,
                           node3}; // store selected node numbers used in n-opt
                                   // int opted_path[7][n];  // store total node
                                   // list of each roots. ith global array.
  double min_path_length;
  int min_path_num = -1; // number of index which gives minimum length
  // -1 means no change
  int i;
  if (whether_geograph == 0) {
    min_path_length =
        euclid_distance(cities, new_path[node1], new_path[node1 + 1]) +
        euclid_distance(cities, new_path[node2], new_path[node2 + 1]) +
        euclid_distance(cities, new_path[node3], new_path[node3 + 1]);
  } else {
    min_path_length =
        geograph_distance(cities, new_path[node1], new_path[node1 + 1]) +
        geograph_distance(cities, new_path[node2], new_path[node2 + 1]) +
        geograph_distance(cities, new_path[node3], new_path[node3 + 1]);
  }

  opt_3_2_cal_length(cities, new_path, n, selected_nodes, path_length,
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
    if (min_path_num == 0) {
      reverse(new_path, n, node1, node3);
      reverse(new_path, n, node2, node3);
    } else if (min_path_num == 1) {
      reverse(new_path, n, node1, node2);
      reverse(new_path, n, node2, node3);
    } else if (min_path_num == 2) {
      reverse(new_path, n, node1, node2);
      reverse(new_path, n, node1, node3);
    } else if (min_path_num == 3) {
      reverse(new_path, n, node1, node2);
    } else if (min_path_num == 4) {
      reverse(new_path, n, node2, node3);
    } else if (min_path_num == 5) {
      reverse(new_path, n, node1, node3);
    } else if (min_path_num == 6) {
      reverse(new_path, n, node1, node3);
      reverse(new_path, n, node2, node3);
    } else if (min_path_num == 3) {
      reverse(new_path, n, node1, node2);
    } else if (min_path_num == 4) {
      reverse(new_path, n, node2, node3);
    } else if (min_path_num == 5) {
      reverse(new_path, n, node1, node3);
    } else if (min_path_num == 6) {
      reverse(new_path, n, node1, node3);
      reverse(new_path, n, node2, node3);
    }
    *new_path_length = distance(cities, new_path, n, whether_geograph);
    return 1;
  } else {
    return 0;
  }
}
int main(void) {
  int N, i, j;
  struct city cities[maxN]; // city's data structure :array
  int new_path[maxN];
  int visited[maxN];
  int path_size = 0;
  int improved = 0; // identify whether toor is improved by n-opt or not;
  // 0: no improve; 1: improved;
  double new_path_length;
  double best_path_length = inf;
  int best_path[maxN];
  int whether_geograph;         // 0: euclidean 1:geographical
  double time_limit = 1.034654; // 計算時間制限（秒）
  clock_t start_t, end_t;
  double utime;
  char fname[128];
  FILE *fp;
  int buf;

  srand(time(NULL)); // initialize rand

  printf("input filename: ");
  scanf("%s", fname);
  fp = fopen(fname, "r");
  printf("input data's type 0:euclidean 1:geographical\n");
  scanf("%d", &whether_geograph);

  if (fp == NULL) {
    printf("Error: Unable to open file.\n");
    return 1;
  }

  fscanf(fp, "%d", &N); // 頂点数を読み込む
  for (i = 0; i < N; i++) {
    fscanf(fp, "%d %lf %lf", &buf, &cities[i].x,
           &cities[i].y); // 座標を読み込む
  }
  fclose(fp);

  // 距離行列の計算
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      if (whether_geograph == 0) {
        dist[i][j] = euclid_distance(cities, i, j);
      } else {
        dist[i][j] = geograph_distance(cities, i, j);
      }
    }
  }

  // 初期化
  for (i = 0; i < N; i++) {
    visited[i] = -1;
  }

  start_t = clock();

  while (1) {
    // 初期解の生成
    new_path[0] = 0; // 任意の初期点（ここでは0を選択）
    visited[0] = 1;
    path_size = 1;

    while (path_size < N) {
      int farthest = find_farthest_node(N, visited);
      insert_nodes(new_path, &path_size, farthest);
      visited[farthest] = 1;
    }

    // 巡回路の総距離の計算
    new_path_length = distance(cities, new_path, N, path_size);

    // 3+2-optによる改善
    for (i = 0; i < 1000; i++) { // 1000回繰り返して改善（適宜調整可能）
      int node1 = createrand(N);
      int node2 = createrand(N);
      int node3 = createrand(N);
      improved = opt3_2_once(cities, new_path, path_size, &new_path_length,
                             node1, node2, node3, whether_geograph);
    }

    end_t = clock();
    utime = (double)(end_t - start_t) / CLOCKS_PER_SEC;

    if (new_path_length < best_path_length) {
      best_path_length = new_path_length;
      for (i = 0; i < path_size; i++) {
        best_path[i] = new_path[i];
      }
    }

    if (utime > time_limit) {
      break;
    }
  }

  // 結果の出力
  printf("Best Path:\n");
  for (i = 0; i < path_size; i++) {
    printf("%d ", best_path[i]);
  }
  printf("\nTotal Distance: %f\n", best_path_length);
  printf("Calculation Time: %f seconds\n", utime);

  return 0;
}
