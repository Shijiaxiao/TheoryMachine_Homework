/* 80 150 200 140 75 25 100 0 1 */

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <opencv2/opencv.hpp>

using std::vector;
using std::string;

const double PI = 3.1415926;

// Define serial index:
const int point_nums = 6;
enum serial { A, B, C, D, E, F };
enum subscript { t, ang, x, y };

// Cordination of each point:
vector<vector<double>> pos_temp(point_nums, vector<double>(4));
vector<vector<vector<double>>> pos(point_nums, vector<vector<double>>(0, vector<double>(4)));

// Linear speed of each point:
vector<vector<double>> v_temp(point_nums, vector<double>(4));
vector<vector<vector<double>>> v(point_nums, vector<vector<double>>(0, vector<double>(4)));

// Acceleration of each point:
vector<vector<double>> a_temp(point_nums, vector<double>(4));
vector<vector<vector<double>>> a(point_nums, vector<vector<double>>(0, vector<double>(4)));

// Angular speed and angle of driving member:
double omega;       // 100
double angle_init;  // 0
double theta_init;  // 0
double angle;
double theta;
double d_angle;
double dt;          // PI / 18000 (s)
double T = 0;

// Length of each pole (mm):
double lab;   // 80
double lbc;   // 150
double lcd;   // 200
double lad;   // 140
double lbe;   // 75
double lef;   // 25

// Transform of double to string:
string formatDoubleValue(double val, int fixed);

// Initialization:
void initialize();

// Calculation:
void calPosPoint(double T, int index);
void calVPoint(serial alpha, double T, int index);
void calAPoint(serial alpha, double T, int index);

// Generate pictures:
void drawPoint(serial alpha);
void drawVLine(serial alpha);
void drawALine(serial alpha);

// Generate movements:
void drawMovements(serial alpha);

int main(int argc, const char * argv[]) {
    
    initialize();
    
    // Calculate position:
    T = 0;
    for (int i = 0; i <= int(360 / d_angle); i++) {
        T = dt * i;
        calPosPoint(T, i);
        // printf("pos_F (%lf, %.1lf, %.3lf, %.3lf)\n", pos[F][i][t], pos[F][i][ang], pos[F][i][x], pos[F][i][y]);
    }
    
    // Calculate vel:
    T = 0;
    for (int i = 0; i <= int(360 / d_angle); i++) {
        T = dt * i;
        calVPoint(F, T, i);
        // printf("vel_F (%lf, %.1lf, %.1lf, %.1lf)\n", v[F][i][t], v[F][i][ang], v[F][i][x], v[F][i][y]);
    }
    
    // Calculate acceleration:
    T = 0;
    for (int i = 0; i <= int(360 / d_angle); i++) {
        T = dt * i;
        calAPoint(F, T, i);
        // printf("acc_F (%lf, %.1lf, %.1lf, %.1lf)\n", a[F][i][t], a[F][i][ang], a[F][i][x], a[F][i][y]);
    }
    
    // Generate pictures:
    drawPoint(F);
    drawVLine(F);
    drawALine(F);
    
    // Generate movements:
    drawMovements(F);
    
    return 0;
}

string formatDoubleValue(double val, int fixed) {
    auto str = std::to_string(val);
    return str.substr(0, str.find(".") + fixed + 1);
}

void initialize() { 
    // Input:
    printf("Please input lab: ");
    scanf("%lf", &lab);
    printf("Please input lbc: ");
    scanf("%lf", &lbc);
    printf("Please input lcd: ");
    scanf("%lf", &lcd);
    printf("Please input lad: ");
    scanf("%lf", &lad);
    printf("Please input lbe: ");
    scanf("%lf", &lbe);
    printf("Please input lef: ");
    scanf("%lf", &lef);
    printf("Please input omega: ");
    scanf("%lf", &omega);
    printf("Please input initial theta_init: ");
    scanf("%lf", &angle_init);
    theta_init = angle_init * PI / 180;
    printf("Please input initial d_angle: ");
    scanf("%lf", &d_angle);
    
    dt = (d_angle * PI) / (omega * 180);
}

void calPosPoint(double T, int index) {
    
    theta = theta_init + omega * T;
    angle = angle_init + index + d_angle; 
    // angle = theta * 180 / PI;
    
    // Point A:
    pos_temp[A][t] = T;
    pos_temp[A][ang] = angle;
    pos_temp[A][x] = 0.0;
    pos_temp[A][y] = 0.0;
    pos[A].push_back(pos_temp[A]);
    
    // Point B:
    pos_temp[B][t] = T;
    pos_temp[B][ang] = angle;
    pos_temp[B][x] = lab * cos(theta);
    pos_temp[B][y] = lab * sin(theta);
    pos[B].push_back(pos_temp[B]);
    
    // Point C:
    double lbd = sqrt(pow(lab, 2) + pow(lad, 2) - 2 * lab * lad * cos(theta));                // Length of BD
    double theta_adb = asin(pos_temp[B][y] / lbd);                                            // Angle ADB
    double theta_bdc = acos((pow(lbd, 2) + pow(lcd, 2) - pow(lbc, 2)) / (2 * lbd * lcd));     // Angle BDC
    double theta_adc = theta_adb + theta_bdc;                                                 // Angle ADC
    pos_temp[C][t] = T;
    pos_temp[C][ang] = angle;
    pos_temp[C][x] = lad - lcd * cos(theta_adc);
    pos_temp[C][y] = lcd * sin(theta_adc);
    pos[C].push_back(pos_temp[C]);
    
    // Point D:
    pos_temp[D][t] = T;
    pos_temp[D][ang] = angle;
    pos_temp[D][x] = lad;
    pos_temp[D][y] = 0.0;
    pos[D].push_back(pos_temp[D]);
    
    // Point E:
    pos_temp[E][t] = T;
    pos_temp[E][ang] = angle;
    pos_temp[E][x] = pos_temp[B][x] + (pos_temp[C][x] - pos_temp[B][x]) * lbe / lbc;
    pos_temp[E][y] = pos_temp[B][y] + (pos_temp[C][y] - pos_temp[B][y]) * lbe / lbc;
    pos[E].push_back(pos_temp[E]);
    
    // Point F:
    double theta_cbx = atan2(pos_temp[C][y] - pos_temp[B][y], pos_temp[C][x] - pos_temp[B][x]);
    double theta_fex = PI / 2 + theta_cbx;
    pos_temp[F][t] = T;
    pos_temp[F][ang] = angle;
    pos_temp[F][x] = pos_temp[E][x] + lef * cos(theta_fex);;
    pos_temp[F][y] = pos_temp[E][y] + lef * sin(theta_fex);
    pos[F].push_back(pos_temp[F]);
    
    // Print position of each point:
    // printf("\nAngle: %lf\n", angle);
    // printf("A (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[A][t], pos_temp[A][ang], pos_temp[A][x], pos_temp[A][y]);
    // printf("B (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[B][t], pos_temp[B][ang], pos_temp[B][x], pos_temp[B][y]);
    // printf("C (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[C][t], pos_temp[C][ang], pos_temp[C][x], pos_temp[C][y]);
    // printf("D (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[D][t], pos_temp[D][ang], pos_temp[D][x], pos_temp[D][y]);
    // printf("E (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[E][t], pos_temp[E][ang], pos_temp[E][x], pos_temp[E][y]);
    // printf("F (%lf, %.1lf, %.1lf, %.1lf)\n", pos_temp[F][t], pos_temp[F][ang], pos_temp[F][x], pos_temp[F][y]);
    // putchar('\n');
}

void calVPoint(serial alpha, double T, int index) {
    
    angle = angle_init + index * d_angle;
    
    double this_x = pos[alpha][index][x];
    double this_y = pos[alpha][index][y];
    double last_x;
    double last_y;

    if (index == 0) {
        last_x = pos[alpha][int(360/d_angle)-1][x];
        last_y = pos[alpha][int(360/d_angle)-1][y];
    } else {
        last_x = pos[alpha][index-1][x];
        last_y = pos[alpha][index-1][y];
    }

    v_temp[alpha][t] = T;
    v_temp[alpha][ang] = angle;
    v_temp[alpha][x] = (this_x - last_x) / dt;
    v_temp[alpha][y] = (this_y - last_y) / dt;
    v[alpha].push_back(v_temp[alpha]);
    
}

void calAPoint(serial alpha, double T, int index) {
    
    angle = angle_init + index * d_angle;
    
    double this_vx = v[alpha][index][x];
    double this_vy = v[alpha][index][y];
    double last_vx;
    double last_vy;

    if (index == 0) {
        last_vx = v[alpha][int(360/d_angle)-1][x];
        last_vy = v[alpha][int(360/d_angle)-1][y];
    } else {
        last_vx = v[alpha][index-1][x];
        last_vy = v[alpha][index-1][y];
    }

    a_temp[alpha][t] = T;
    a_temp[alpha][ang] = angle;
    a_temp[alpha][x] = (this_vx - last_vx) / dt;
    a_temp[alpha][y] = (this_vy - last_vy) / dt;
    a[alpha].push_back(a_temp[alpha]);

}

void drawPoint(serial alpha) {
    
    vector<cv::Point2f> position;

    // Procession of x and y:
    const int x_count = 8;
    const int y_count = 8;
    const int x_add = 150;
    const int y_add = 150;

    double x_max = 0;
    int x_max_index = 0;

    double x_min = 0;
    int x_min_index = 0;

    double y_max = 0;
    int y_max_index = 0;

    double y_min = 0;
    int y_min_index = 0;

    // Put all points into vector:
    for (int i = 0; i < pos[alpha].size(); i++) {
        position.push_back(cv::Point2f(pos[alpha][i][x] * x_count, pos[alpha][i][y] * y_count));
        if (pos[alpha][i][x] > x_max) { x_max = pos[alpha][i][x]; x_max_index = i; }
        if (pos[alpha][i][x] < x_min) { x_min = pos[alpha][i][x]; x_min_index = i; }
        if (pos[alpha][i][y] > y_max) { y_max = pos[alpha][i][y]; y_max_index = i; }
        if (pos[alpha][i][y] < y_min) { y_min = pos[alpha][i][y]; y_min_index = i; }
    }
    
    // Update max and min:
    x_max *= x_count;
    x_min *= x_count;
    y_max *= y_count;
    y_min *= y_count;

    printf("x_max: %d, x_min: %d, y_max: %d, vy_min: %d\n", int(x_max / x_count), int(x_min / x_count), int(y_max / y_count), int(y_min / y_count));
    printf("x_range: %d, y_range: %d\n", int((x_max - x_min) / x_count), int((y_max - y_min) / y_count));

    // Create image and initialize it: 
    cv::Mat image_position = cv::Mat::zeros(int(y_max-y_min)+4*y_add, int(x_max-x_min)+4*x_add, CV_8UC3);

    // Draw x and y line:
    cv::line(image_position, cv::Point(x_add, int(y_max)+2*y_add), cv::Point(image_position.cols-x_add, int(y_max)+2*y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add, y_add), cv::Point(int(-x_min)+2*x_add, image_position.rows-y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_position, "0", cv::Point(int(-x_min)+2*x_add-55, int(y_max)+2*y_add+55), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_position, "Fx", cv::Point(image_position.cols-x_add+40, int(y_max)+2*y_add+18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_position, "Fy", cv::Point(int(-x_min)+2*x_add-28, y_add-40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_position, cv::Point(int(-x_min)+2*x_add+int(pos.x), int(y_max)+2*y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // x_max
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[x_max_index].x, int(y_max)+2*y_add-position[x_max_index].y), cv::Point(int(-x_min)+2*x_add, int(y_max)+2*y_add-position[x_max_index].y), cv::Scalar(230,216,173), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[x_max_index].x, int(y_max)+2*y_add-position[x_max_index].y), cv::Point(int(-x_min)+2*x_add+position[x_max_index].x, int(y_max)+2*y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_position, formatDoubleValue(x_max/x_count, 1), cv::Point(int(-x_min)+2*x_add+int(x_max)-50, int(y_max)+2*y_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[x_max_index].y/y_count, 1), cv::Point(int(-x_min)+2*x_add-150, int(y_max)+2*y_add-position[x_max_index].y+70), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // x_min
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[x_min_index].x, int(y_max)+2*y_add-position[x_min_index].y), cv::Point(int(-x_min)+2*x_add, int(y_max)+2*y_add-position[x_min_index].y), cv::Scalar(230,216,173), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[x_min_index].x, int(y_max)+2*y_add-position[x_min_index].y), cv::Point(int(-x_min)+2*x_add+position[x_min_index].x, int(y_max)+2*y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_position, formatDoubleValue(x_min/x_count, 1), cv::Point(int(-x_min)+2*x_add+int(x_min)-80, int(y_max)+2*y_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[x_min_index].y/y_count, 1), cv::Point(int(-x_min)+2*x_add+15, int(y_max)+2*y_add-position[x_min_index].y-10), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    
    // y_max
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[y_max_index].x, int(y_max)+2*y_add-position[y_max_index].y), cv::Point(int(-x_min)+2*x_add, int(y_max)+2*y_add-position[y_max_index].y), cv::Scalar(230,216,173), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[y_max_index].x, int(y_max)+2*y_add-position[y_max_index].y), cv::Point(int(-x_min)+2*x_add+position[y_max_index].x, int(y_max)+2*y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_position, formatDoubleValue(y_max/y_count, 1), cv::Point(int(-x_min)+2*x_add-150, 2*y_add+12), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[y_max_index].x/x_count, 1), cv::Point(int(-x_min)+2*x_add+position[y_max_index].x-50, int(y_max)+2*y_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    
    // y_min
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[y_min_index].x, int(y_max)+2*y_add-position[y_min_index].y), cv::Point(int(-x_min)+2*x_add, int(y_max)+2*y_add-position[y_min_index].y), cv::Scalar(230,216,173), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[y_min_index].x, int(y_max)+2*y_add-position[y_min_index].y), cv::Point(int(-x_min)+2*x_add+position[y_min_index].x, int(y_max)+2*y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_position, formatDoubleValue(y_min/y_count, 1), cv::Point(int(-x_min)+2*x_add+15, int(y_max-y_min)+2*y_add+12), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[y_min_index].x/x_count, 1), cv::Point(int(-x_min)+2*x_add+position[y_min_index].x-70, int(y_max)+2*y_add-28), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    
    // Draw start point:
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[0].x, int(y_max)+2*y_add-position[0].y), cv::Point(int(-x_min)+2*x_add+position[0].x, int(y_max)+2*y_add), cv::Scalar(160,160,255), 1);
    cv::line(image_position, cv::Point(int(-x_min)+2*x_add+position[0].x, int(y_max)+2*y_add-position[0].y), cv::Point(int(-x_min)+2*x_add, int(y_max)+2*y_add-position[0].y), cv::Scalar(160,160,255), 1);
    cv::putText(image_position, "start", cv::Point(int(-x_min)+2*x_add+position[0].x+35, int(y_max)+2*y_add-position[0].y+18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[0].x/x_count, 1), cv::Point(int(-x_min)+2*x_add+position[0].x-30, int(y_max)+2*y_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(160,160,255), 1, 6);
    cv::putText(image_position, formatDoubleValue(position[0].y/y_count, 1), cv::Point(int(-x_min)+2*x_add-120, int(y_max)+2*y_add-position[0].y+20), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(160,160,255), 1, 6);
    
    // Show images:
    // cv::imshow("Image_position", image_position);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/Image_position_F.tiff", image_position);
    cv::imwrite("../images/JPEG/Image_position_F.jpeg", image_position);
    
}

void drawVLine(serial alpha) {

    vector<cv::Point2f> points_vx_angle;
    vector<cv::Point2f> points_vy_angle;

    // Procession of v:
    const int vx_count = 10;
    const int vy_count = 20;
    const int v_add = 150;

    // Procession of angle:
    const int angle_count = 10;
    const int angle_add1 = 100;
    const int angle_add2 = 300;
    const int angle_offset = 220;

    double vx_max = 0;
    int vx_max_index = 0;

    double vx_min = 0;
    int vx_min_index = 0;

    double vy_max = 0;
    int vy_max_index = 0;

    double vy_min = 0;
    int vy_min_index = 0;

    int angle_max = int(v[alpha].size() * angle_count);

    // Put points into vector:
    for (int i = 0; i < v[alpha].size(); i++) {
        points_vx_angle.push_back(cv::Point2f(v[alpha][i][ang]*angle_count/d_angle+angle_offset, v[alpha][i][x]/vx_count));
        points_vy_angle.push_back(cv::Point2f(v[alpha][i][ang]*angle_count/d_angle+angle_offset, v[alpha][i][y]/vy_count));
        if (v[alpha][i][x] > vx_max) { vx_max = v[alpha][i][x]; vx_max_index = i; }
        if (v[alpha][i][x] < vx_min) { vx_min = v[alpha][i][x]; vx_min_index = i; }
        if (v[alpha][i][y] > vy_max) { vy_max = v[alpha][i][y]; vy_max_index = i; }
        if (v[alpha][i][y] < vy_min) { vy_min = v[alpha][i][y]; vy_min_index = i; }
    }
    
    // Update max and min:
    vx_max /= vx_count;
    vx_min /= vx_count;
    vy_max /= vy_count;
    vy_min /= vy_count;

    printf("vx_max: %d, vx_min: %d, vy_max: %d, vy_min: %d\n", int(vx_max*vx_count/10), int(vx_min*vx_count/10), int(vy_max*vy_count/10), int(vy_min*vy_count/10));
    printf("angle: %d, vx_range: %d, vy_range: %d\n", angle_max, int((vx_max-vx_min)*vx_count/10), int((vy_max-vy_min)*vy_count/10));

    // Create mat to draw function:  angle v(cm/s)
    cv::Mat image_vx_angle = cv::Mat::zeros(int(vx_max-vx_min+4*v_add), angle_max+angle_add1+angle_add2+angle_offset, CV_8UC3);
    cv::Mat image_vy_angle = cv::Mat::zeros(int(vy_max-vy_min+4*v_add), angle_max+angle_add1+angle_add2+angle_offset, CV_8UC3);

    // Draw v angle line:
    cv::line(image_vx_angle, cv::Point(angle_offset, int(vx_max)+2*v_add), cv::Point(image_vx_angle.cols-angle_add2, int(vx_max)+2*v_add), cv::Scalar(255,255,255), 1);
    cv::line(image_vx_angle, cv::Point(angle_offset, int(vx_max-vx_min+3*v_add)), cv::Point(angle_offset, v_add), cv::Scalar(255,255,255), 1);
    cv::line(image_vy_angle, cv::Point(angle_offset, int(vy_max)+2*v_add), cv::Point(image_vy_angle.cols-angle_add2, int(vy_max)+2*v_add), cv::Scalar(255,255,255), 1);
    cv::line(image_vy_angle, cv::Point(angle_offset, int(vy_max-vy_min+3*v_add)), cv::Point(angle_offset, v_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point and texts:
    cv::putText(image_vx_angle, "0", cv::Point(angle_offset-55, int(vx_max)+2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_vx_angle, "360", cv::Point(image_vx_angle.cols-angle_add1-angle_add2-50, int(vx_max)+2*v_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vx_angle, "angle", cv::Point(image_vx_angle.cols-angle_add2+40, int(vx_max)+2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_vx_angle, "Fvx(cm/s)", cv::Point(angle_offset-110, v_add-50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    cv::putText(image_vy_angle, "0", cv::Point(angle_offset-55, int(vy_max)+2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_vy_angle, "360", cv::Point(image_vy_angle.cols-angle_add1-angle_add2-50, int(vy_max)+2*v_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vy_angle, "angle", cv::Point(image_vy_angle.cols-angle_add2+40, int(vy_max)+2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_vy_angle, "Fvy(cm/s)", cv::Point(angle_offset-110, v_add-50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);

    // Draw points on mat:
    for (int i = 0; i < points_vx_angle.size(); i++) {
        cv::Point2f pt_vx_angle = points_vx_angle[i];
        cv::Point2f pt_vy_angle = points_vy_angle[i];
        cv::circle(image_vx_angle, cv::Point(int(pt_vx_angle.x), int(vx_max)+2*v_add - int(pt_vx_angle.y)), 3, cv::Scalar(10,215,255), -1);
        cv::circle(image_vy_angle, cv::Point(int(pt_vy_angle.x), int(vy_max)+2*v_add - int(pt_vy_angle.y)), 3, cv::Scalar(10,215,255), -1);
    }

    // vx_max
    cv::line(image_vx_angle, cv::Point(points_vx_angle[vx_max_index].x, int(vx_max)+2*v_add - int(points_vx_angle[vx_max_index].y)), cv::Point(points_vx_angle[vx_max_index].x, int(vx_max)+2*v_add), cv::Scalar(230,216,173), 1);
    cv::line(image_vx_angle, cv::Point(points_vx_angle[vx_max_index].x, int(vx_max)+2*v_add - int(points_vx_angle[vx_max_index].y)), cv::Point(angle_offset, int(vx_max)+2*v_add - int(points_vx_angle[vx_max_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_vx_angle, std::to_string(int(vx_max*vx_count/10)), cv::Point(angle_offset-115, 2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vx_angle, std::to_string(int(vx_max_index*d_angle)), cv::Point(points_vx_angle[vx_max_index].x-32, int(vx_max)+2*v_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // vx_min
    cv::line(image_vx_angle, cv::Point(points_vx_angle[vx_min_index].x, int(vx_max)+2*v_add - int(points_vx_angle[vx_min_index].y)), cv::Point(points_vx_angle[vx_min_index].x, int(vx_max)+2*v_add), cv::Scalar(230,216,173), 1);
    cv::line(image_vx_angle, cv::Point(points_vx_angle[vx_min_index].x, int(vx_max)+2*v_add - int(points_vx_angle[vx_min_index].y)), cv::Point(angle_offset, int(vx_max)+2*v_add - int(points_vx_angle[vx_min_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_vx_angle, std::to_string(int(vx_min*vx_count/10)), cv::Point(angle_offset-145, int(vx_max-vx_min+2*v_add+15)), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vx_angle, std::to_string(int(vx_min_index*d_angle)), cv::Point(points_vx_angle[vx_min_index].x-50, int(vx_max)+2*v_add-25), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // vx 360
    cv::line(image_vx_angle, cv::Point(points_vx_angle[v[alpha].size()-1].x, int(vx_max)+2*v_add - int(points_vx_angle[v[alpha].size()-1].y)), cv::Point(points_vx_angle[v[alpha].size()-1].x, int(vx_max)+2*v_add), cv::Scalar(230,216,173), 1);

    // vy_max
    cv::line(image_vy_angle, cv::Point(points_vy_angle[vy_max_index].x, int(vy_max)+2*v_add - int(points_vy_angle[vy_max_index].y)), cv::Point(points_vy_angle[vy_max_index].x, int(vy_max)+2*v_add), cv::Scalar(230,216,173), 1);
    cv::line(image_vy_angle, cv::Point(points_vy_angle[vy_max_index].x, int(vy_max)+2*v_add - int(points_vy_angle[vy_max_index].y)), cv::Point(angle_offset, int(vy_max)+2*v_add - int(points_vy_angle[vy_max_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_vy_angle, std::to_string(int(vy_max*vy_count/10)), cv::Point(angle_offset-145, 2*v_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vy_angle, std::to_string(int(vy_max_index*d_angle)), cv::Point(points_vy_angle[vy_max_index].x-16, int(vy_max)+2*v_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // vy_min
    cv::line(image_vy_angle, cv::Point(points_vy_angle[vy_min_index].x, int(vy_max)+2*v_add - int(points_vy_angle[vy_min_index].y)), cv::Point(points_vy_angle[vy_min_index].x, int(vy_max)+2*v_add), cv::Scalar(230,216,173), 1);
    cv::line(image_vy_angle, cv::Point(points_vy_angle[vy_min_index].x, int(vy_max)+2*v_add - int(points_vy_angle[vy_min_index].y)), cv::Point(angle_offset, int(vy_max)+2*v_add - int(points_vy_angle[vy_min_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_vy_angle, std::to_string(int(vy_min*vy_count/10)), cv::Point(angle_offset-145, int(vy_max-vy_min+2*v_add+15)), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_vy_angle, std::to_string(int(vy_min_index*d_angle)), cv::Point(points_vy_angle[vy_min_index].x-50, int(vy_max)+2*v_add-25), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // vy 360
    cv::line(image_vy_angle, cv::Point(points_vy_angle[v[alpha].size()-1].x, int(vy_max)+2*v_add - int(points_vy_angle[v[alpha].size()-1].y)), cv::Point(points_vy_angle[v[alpha].size()-1].x, int(vy_max)+2*v_add), cv::Scalar(230,216,173), 1);

    // Show images:
    // cv::imshow("Image_vx_angle", image_vx_angle);
    // cv::imshow("Image_vy_angle", image_vy_angle);
    // cv::waitKey();

    // Save images:
    cv::imwrite("../images/TIFF/Image_Fvx_angle.tiff", image_vx_angle);
    cv::imwrite("../images/TIFF/Image_Fvy_angle.tiff", image_vy_angle);
    cv::imwrite("../images/JPEG/Image_Fvx_angle.jpeg", image_vx_angle);
    cv::imwrite("../images/JPEG/Image_Fvy_angle.jpeg", image_vy_angle);
    
}

void drawALine(serial alpha) {

    vector<cv::Point2f> points_ax_angle;
    vector<cv::Point2f> points_ay_angle;

    // Procession of a:
    const int ax_count = 3000;
    const int ay_count = 4000;
    const int a_add = 150;

    // Procession of angle:
    const int angle_count = 10;
    const int angle_add1 = 100;
    const int angle_add2 = 300;
    const int angle_offset = 240;

    double ax_max = 0;
    int ax_max_index = 0;

    double ax_min = 0;
    int ax_min_index = 0;

    double ay_max = 0;
    int ay_max_index = 0;

    double ay_min = 0;
    int ay_min_index = 0;

    int angle_max = int(a[alpha].size() * angle_count);

    // Put all points into vector:
    for (int i = 0; i < a[alpha].size(); i++) {
        points_ax_angle.push_back(cv::Point2f(a[alpha][i][ang]*angle_count/d_angle+angle_offset, a[alpha][i][x]/ax_count));
        points_ay_angle.push_back(cv::Point2f(a[alpha][i][ang]*angle_count/d_angle+angle_offset, a[alpha][i][y]/ay_count));
        if (a[alpha][i][x] > ax_max) { ax_max = a[alpha][i][x]; ax_max_index = i; }
        if (a[alpha][i][x] < ax_min) { ax_min = a[alpha][i][x]; ax_min_index = i; }
        if (a[alpha][i][y] > ay_max) { ay_max = a[alpha][i][y]; ay_max_index = i; }
        if (a[alpha][i][y] < ay_min) { ay_min = a[alpha][i][y]; ay_min_index = i; }
    }
    
    // Update max and min:
    ax_max /= ax_count;
    ax_min /= ax_count;
    ay_max /= ay_count;
    ay_min /= ay_count;

    printf("ax_max: %d, ax_min: %d, ay_max: %d, ay_min: %d\n", int(ax_max*ax_count/1000), int(ax_min*ax_count/1000), int(ay_max*ay_count/1000), int(ay_min*ay_count/1000));
    printf("angle: %d, ax_range: %d, ay_range: %d\n", angle_max, int((ax_max-ax_min)*ax_count/1000), int((ay_max-ay_min)*ay_count/1000));

    // Create mat to draw function:  angle v: cm/s
    cv::Mat image_ax_angle = cv::Mat::zeros(int(ax_max-ax_min+4*a_add), angle_max+angle_add1+angle_add2+angle_offset, CV_8UC3);
    cv::Mat image_ay_angle = cv::Mat::zeros(int(ay_max-ay_min+4*a_add), angle_max+angle_add1+angle_add2+angle_offset, CV_8UC3);

    // Draw v angle line:
    cv::line(image_ax_angle, cv::Point(angle_offset, int(ax_max)+2*a_add), cv::Point(image_ax_angle.cols-angle_add2, int(ax_max)+2*a_add), cv::Scalar(255,255,255), 1);
    cv::line(image_ax_angle, cv::Point(angle_offset, int(ax_max-ax_min+3*a_add)), cv::Point(angle_offset, a_add), cv::Scalar(255,255,255), 1);
    cv::line(image_ay_angle, cv::Point(angle_offset, int(ay_max)+2*a_add), cv::Point(image_ay_angle.cols-angle_add2, int(ay_max)+2*a_add), cv::Scalar(255,255,255), 1);
    cv::line(image_ay_angle, cv::Point(angle_offset, int(ay_max-ay_min+3*a_add)), cv::Point(angle_offset, a_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point and texts:
    cv::putText(image_ax_angle, "0", cv::Point(angle_offset-55, int(ax_max)+2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_ax_angle, "360", cv::Point(image_ax_angle.cols-angle_add1-angle_add2-50, int(ax_max)+2*a_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ax_angle, "angle", cv::Point(image_ax_angle.cols-angle_add2+40, int(ax_max)+2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_ax_angle, "Fax(m/s2)", cv::Point(angle_offset-110, a_add-50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    cv::putText(image_ay_angle, "0", cv::Point(angle_offset-55, int(ay_max)+2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_ay_angle, "360", cv::Point(image_ay_angle.cols-angle_add1-angle_add2-20, int(ay_max)+2*a_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ay_angle, "angle", cv::Point(image_ay_angle.cols-angle_add2+40, int(ay_max)+2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_ay_angle, "Fay(m/s2)", cv::Point(angle_offset-110, a_add-50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);

    // Draw points on mat:
    for (int i = 0; i < points_ax_angle.size(); i++) {
        cv::Point2f pt_vx_angle = points_ax_angle[i];
        cv::Point2f pt_vy_angle = points_ay_angle[i];
        cv::circle(image_ax_angle, cv::Point(int(pt_vx_angle.x), int(ax_max)+2*a_add - int(pt_vx_angle.y)), 3, cv::Scalar(10,215,255), -1);
        cv::circle(image_ay_angle, cv::Point(int(pt_vy_angle.x), int(ay_max)+2*a_add - int(pt_vy_angle.y)), 3, cv::Scalar(10,215,255), -1);
    }

    // ax_max
    cv::line(image_ax_angle, cv::Point(points_ax_angle[ax_max_index].x, int(ax_max)+2*a_add - int(points_ax_angle[ax_max_index].y)), cv::Point(points_ax_angle[ax_max_index].x, int(ax_max)+2*a_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ax_angle, cv::Point(points_ax_angle[ax_max_index].x, int(ax_max)+2*a_add - int(points_ax_angle[ax_max_index].y)), cv::Point(angle_offset, int(ax_max)+2*a_add - int(points_ax_angle[ax_max_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_ax_angle, std::to_string(int(ax_max*ax_count/1000)), cv::Point(angle_offset-155, 2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ax_angle, std::to_string(int(ax_max_index*d_angle)), cv::Point(points_ax_angle[ax_max_index].x-20, int(ax_max)+2*a_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // ax_min
    cv::line(image_ax_angle, cv::Point(points_ax_angle[ax_min_index].x, int(ax_max)+2*a_add - int(points_ax_angle[ax_min_index].y)), cv::Point(points_ax_angle[ax_min_index].x, int(ax_max)+2*a_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ax_angle, cv::Point(points_ax_angle[ax_min_index].x, int(ax_max)+2*a_add - int(points_ax_angle[ax_min_index].y)), cv::Point(angle_offset, int(ax_max)+2*a_add - int(points_ax_angle[ax_min_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_ax_angle, std::to_string(int(ax_min*ax_count/1000)), cv::Point(angle_offset-185, int(ax_max-ax_min+2*a_add+15)), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ax_angle, std::to_string(int(ax_min_index*d_angle)), cv::Point(points_ax_angle[ax_min_index].x-30, int(ax_max)+2*a_add-25), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // ax 360 
    cv::line(image_ax_angle, cv::Point(points_ax_angle[v[alpha].size()-1].x, int(ax_max)+2*a_add - int(points_ax_angle[v[alpha].size()-1].y)), cv::Point(points_ax_angle[v[alpha].size()-1].x, int(ax_max)+2*a_add), cv::Scalar(230,216,173), 1);

    // ay_max
    cv::line(image_ay_angle, cv::Point(points_ay_angle[ay_max_index].x, int(ay_max)+2*a_add - int(points_ay_angle[ay_max_index].y)), cv::Point(points_ay_angle[ay_max_index].x, int(ay_max)+2*a_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ay_angle, cv::Point(points_ay_angle[ay_max_index].x, int(ay_max)+2*a_add - int(points_ay_angle[ay_max_index].y)), cv::Point(angle_offset, int(ay_max)+2*a_add - int(points_ay_angle[ay_max_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_ay_angle, std::to_string(int(ay_max*ay_count/1000)), cv::Point(angle_offset-155, 2*a_add+15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ay_angle, std::to_string(int(ay_max_index*d_angle)), cv::Point(points_ay_angle[ay_max_index].x-80, int(ay_max)+2*a_add+60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // ay_min
    cv::line(image_ay_angle, cv::Point(points_ay_angle[ay_min_index].x, int(ay_max)+2*a_add - int(points_ay_angle[ay_min_index].y)), cv::Point(points_ay_angle[ay_min_index].x, int(ay_max)+2*a_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ay_angle, cv::Point(points_ay_angle[ay_min_index].x, int(ay_max)+2*a_add - int(points_ay_angle[ay_min_index].y)), cv::Point(angle_offset, int(ay_max)+2*a_add - int(points_ay_angle[ay_min_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_ay_angle, std::to_string(int(ay_min*ay_count/1000)), cv::Point(angle_offset-185, int(ay_max-ay_min+2*a_add+15)), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ay_angle, std::to_string(int(ay_min_index*d_angle)), cv::Point(points_ay_angle[ay_min_index].x-30, int(ay_max)+2*a_add-25), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // ay 360
    cv::line(image_ay_angle, cv::Point(points_ay_angle[v[alpha].size()-1].x, int(ay_max)+2*a_add - int(points_ay_angle[v[alpha].size()-1].y)), cv::Point(points_ay_angle[v[alpha].size()-1].x, int(ay_max)+2*a_add), cv::Scalar(230,216,173), 1);

    // Show images:
    // cv::imshow("Image_ax_angle", image_ax_angle);
    // cv::imshow("Image_ay_angle", image_ay_angle);
    // cv::waitKey();

    // Save images:
    cv::imwrite("../images/TIFF/Image_Fax_angle.tiff", image_ax_angle);
    cv::imwrite("../images/TIFF/Image_Fay_angle.tiff", image_ay_angle);
    cv::imwrite("../images/JPEG/Image_Fax_angle.jpeg", image_ax_angle);
    cv::imwrite("../images/JPEG/Image_Fay_angle.jpeg", image_ay_angle);
    
}

void drawMovements(serial alpha) {

    // Initialize MP4V codec
    int fourcc = cv::VideoWriter::fourcc('m', 'p', '4', 'v'); 
    double fps = 20.0;

    // Create images:
    vector<cv::Mat> images;

    // Procession of x and y:
    const int x_count = 6;
    const int y_count = 6;
    const int x_add = 250;
    const int y_add = 250;

    double x_max = 0;
    int x_max_index = 0;

    double x_min = 0;
    int x_min_index = 0;

    double y_max = 0;
    int y_max_index = 0;

    double y_min = 0;
    int y_min_index = 0;

    // Put all points into vector:
    for (int i = 0; i < pos[alpha].size(); i++) {
        for (int j = A; j < F; j++) {
            // printf("\n%d\n", debug++);
            if (pos[j][i][x] > x_max) {
                x_max = pos[j][i][x];
                x_max_index = i * 10 + j;
            }
            if (pos[j][i][x] < x_min) {
                x_min = pos[j][i][x];
                x_min_index = i * 10 + j;
            }
            if (pos[j][i][y] > y_max) {
                y_max = pos[j][i][y];
                y_max_index = i * 10 + j;
            }
            if (pos[j][i][y] < y_min) {
                y_min = pos[j][i][y];
                y_min_index = i * 10 + j;
            }
        }
    }
    
    // Update max and min:
    x_max *= x_count;
    x_min *= x_count;
    y_max *= y_count;
    y_min *= y_count;

    printf("x_max: %d, x_min: %d, y_max: %d, vy_min: %d\n", int(x_max / x_count), int(x_min / x_count), int(y_max / y_count), int(y_min / y_count));
    printf("x_range: %d, y_range: %d\n", int((x_max - x_min) / x_count), int((y_max - y_min) / y_count));

    // Create Vedio Writer:
    cv::Mat image_temp = cv::Mat::zeros(int(y_max-y_min)+4*y_add, int(x_max-x_min)+4*x_add, CV_8UC3);
    cv::VideoWriter video("../images/VIDEO/Movement.mp4", fourcc, fps, image_temp.size());
    
    // Generate video:
    for (int i = 0; i < pos[alpha].size(); i++) {

        // Create image:
        cv::Mat image = cv::Mat::zeros(int(y_max-y_min)+4*y_add, int(x_max-x_min)+4*x_add, CV_8UC3);

        // Put name and author
        cv::putText(image, "Movement", cv::Point(image.cols-2*x_add-200, y_add-30), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(230,216,173), 3, 6);
        cv::putText(image, "Author: ShiJiaxiao", cv::Point(image.cols-2*x_add-320, y_add+100), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(230,216,173), 3, 6);

        {   // Draw initial point(0) and texts:
            cv::line(image, cv::Point(x_add, int(y_max)+2*y_add), cv::Point(image.cols-x_add, int(y_max)+2*y_add), cv::Scalar(255,255,255), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add, y_add), cv::Point(int(-x_min)+2*x_add, image.rows-y_add), cv::Scalar(255,255,255), 2);
            cv::putText(image, "0", cv::Point(int(-x_min)+2*x_add-80, int(y_max)+2*y_add+80), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(255,255,255), 2, 6);
            cv::putText(image, "x", cv::Point(image.cols-x_add+40, int(y_max)+2*y_add+20), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
            cv::putText(image, "y", cv::Point(int(-x_min)+2*x_add-20, y_add-60), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
        }
        {   // Draw machinery:
            
            // Lab Lbc Lcd Lef:
            cv::line(image, cv::Point(int(-x_min)+2*x_add+pos[A][i][x]*x_count, int(y_max)+2*y_add-pos[A][i][y]*y_count), cv::Point(int(-x_min)+2*x_add+pos[B][i][x]*x_count, int(y_max)+2*y_add-pos[B][i][y]*y_count), cv::Scalar(238,130,238), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add+pos[B][i][x]*x_count, int(y_max)+2*y_add-pos[B][i][y]*y_count), cv::Point(int(-x_min)+2*x_add+pos[C][i][x]*x_count, int(y_max)+2*y_add-pos[C][i][y]*y_count), cv::Scalar(238,130,238), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add+pos[C][i][x]*x_count, int(y_max)+2*y_add-pos[C][i][y]*y_count), cv::Point(int(-x_min)+2*x_add+pos[D][i][x]*x_count, int(y_max)+2*y_add-pos[D][i][y]*y_count), cv::Scalar(238,130,238), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add+pos[E][i][x]*x_count, int(y_max)+2*y_add-pos[E][i][y]*y_count), cv::Point(int(-x_min)+2*x_add+pos[F][i][x]*x_count, int(y_max)+2*y_add-pos[F][i][y]*y_count), cv::Scalar(238,130,238), 2);
            // Point A B C D E F
            cv::putText(image, "A", cv::Point(int(-x_min)+2*x_add+pos[A][i][x]*x_count-70, int(y_max)+2*y_add-pos[A][i][y]*y_count-40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
            cv::putText(image, "B", cv::Point(int(-x_min)+2*x_add+pos[B][i][x]*x_count-10, int(y_max)+2*y_add-pos[B][i][y]*y_count-40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
            cv::putText(image, "C", cv::Point(int(-x_min)+2*x_add+pos[C][i][x]*x_count-50, int(y_max)+2*y_add-pos[C][i][y]*y_count-40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
            cv::putText(image, "D", cv::Point(int(-x_min)+2*x_add+pos[D][i][x]*x_count+5 , int(y_max)+2*y_add-pos[D][i][y]*y_count-40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
            cv::putText(image, "E", cv::Point(int(-x_min)+2*x_add+pos[E][i][x]*x_count-50, int(y_max)+2*y_add-pos[E][i][y]*y_count-50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
            cv::putText(image, "F", cv::Point(int(-x_min)+2*x_add+pos[F][i][x]*x_count-125, int(y_max)+2*y_add-pos[F][i][y]*y_count-30), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(230,180,173), 3, 6);

        }

        // Draw track of point alpha
        for (int j = 0; j < i; j++) 
            cv::circle(image, cv::Point(int(-x_min)+2*x_add+pos[alpha][j][x]*x_count, int(y_max)+2*y_add-pos[alpha][j][y]*y_count), 3, cv::Scalar(10,215,255), -1);
        

        // Add image to video:
        images.push_back(image);
        video.write(image);

    }

    for (int i = 0 ; i < 100; i++) {
        video.write(images[images.size()-1]);
    }

    // Release VideoWriter:
    video.release();

}

