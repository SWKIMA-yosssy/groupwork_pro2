/*mainで必要なもの*/
//dist[MaxN][MaxN]重み行列
//visited[MaxN]

/*distとvisitedの初期化*/
//for(int i = 0;i < N ;i++){for(int j;j<N;j++){dist[i][j] == -1}}
//for (int i = 0; i < N; i++){ visited[i] = -1}

/*呼び出し*/
//while (path_size){int farthest = find_farthest_node(dist,N<visited);}
//insert_node(new_path,&path_size,dist,farthest);
//visited[farthest]=1;}


// 最も遠い点を見つける関数
int find_farthest_node(double dist[MAXN][MAXN], int N, int visited[MAXN])
{
    int farthest = -1;
    double max_dist = -1;
    for (int i = 0; i < N; i++)
    {
        if (visited[i] == -1){
            for (int  j = 0; j < N; j++)
            {
                if (visited[j] == -1 && dist[i][j] > max_dist){
                    max_dist = dist[i][j];
                    farthest = i;
                }
            }
        }
    }
    return farthest;
}

// 点をツアーに挿入する関数
void insert_nodes(int new_path[MAXN], int *path_size, double dist[MAXN][MAXN], int node) {
    int best_position = -1;
    double best_cost = DBL_MAX;

    for (int i = 0; i < *path_size ; i++) {
        int j = (i + 1) % *path_size;
        double cost = dist[new_path[i]][node] + dist[node][new_path[j]] - dist[new_path[i]][new_path[j]];
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
