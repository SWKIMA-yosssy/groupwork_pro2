    int N, i, j;
    int new_path[MAXN];
    int path_size = 0;
    double new_path_length;
    double best_path_length = INF;
    int best_path[MAXN];
    double time_limit = 60.0; // 計算時間制限（秒）
    clock_t start_t, end_t;
    double utime;
    char fname[128];
    FILE *fp;

    printf("input filename: ");
    scanf("%s", fname);
    fp = fopen(fname, "r");

    if (fp == NULL) {
        printf("Error: Unable to open file.\n");
        return 1;
    }

    fscanf(fp, "%d", &N); // 頂点数を読み込む
    for (i = 0; i < N; i++) {
        fscanf(fp, "%d %lf %lf", &cities[i].id, &cities[i].x, &cities[i].y); // 座標を読み込む
    }
    fclose(fp);

    // 距離行列の計算
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            dist[i][j] = euclid_distance(i, j);
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
            int farthest = find_farthest_node(dist, N, visited);
            insert_nodes(new_path, &path_size, dist, farthest);
            visited[farthest] = 1;
        }

        // 巡回路の総距離の計算
        new_path_length = distance(new_path, path_size);

        // 3+2-optによる改善
        for (i = 0; i < 1000; i++) { // 1000回繰り返して改善（適宜調整可能）
            int node1 = createrand(N);
            int node2 = createrand(N);
            int node3 = createrand(N);
            opt3_2_once(new_path, path_size, &new_path_length, node1, node2, node3);
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
