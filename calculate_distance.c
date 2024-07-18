#define maxN 500  

/*ユークリッド*/
struct city{
    int x;
    int y;
}

city path[maxN];

double euclid_distance(path[i], path[i+1]) {
    double dx = path[i].x - path[i+1].x;
    double dy = path[i].y - path[i+1].y;
    return sqrt(dx * dx + dy * dy);
}

/*地理的距離*/
#define earthr 6378.388 // 地球の半径

// 度をラジアンに変換する関数
double deg2rad(double deg) {
    return deg * (M_PI / 180.0);
}

/*地理的距離を求める関数*/
double haversine(double lat1, double lon1, double lat2, double lon2) {
    double dlat = deg2rad(lat2 - lat1);
    double dlon = deg2rad(lon2 - lon1);
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(deg2rad(lat1)) * cos(deg2rad(lat2)) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return earthr * c;
}
