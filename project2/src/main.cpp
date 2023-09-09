#include <iostream>
#include <cmath>
#include <vector>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

const double PI = 3.14159265358979323846;
const int    N  = 1000; 

// Designing parameters:
const double stroke = 50.0;              // Stroke
const double cam_angle1 = 130.0;         // Push up
const double cam_angle2 = 100.0;         // Push down
const double pressure_angle1 = 35.0;     // Pressure1
const double pressure_angle2 = 70.0;     // Pressure2
const double stop_angle_far = 80.0;      // Stop far
const double stop_angle_near = 50.0;     // Stop near
const double w = 2 * PI;                 // omega

// Outline parameters:
const double r_roller = 10.0;                          // Roller radius
const double r_base = 70.0;                            // Base radius
const double s0 = 65.0;                                // Base displacement
const double e = sqrt(pow(r_base, 2) - pow(s0, 2));    // e
double min_radius = r_base;                            // Minimum radius

// Movement law:
vector<double> theta(N);
vector<double> theta_pi(N);
vector<double> s(N);
vector<double> v(N);
vector<double> a(N);
vector<double> ds_dtheta(N);
vector<double> d2s_dtheta2(N);
vector<double> alpha(N);
vector<double> radius(N);
const double t1 = cam_angle1 * PI / 180 / (2 * w);
const double t2 = cam_angle2 * PI / 180 / (2 * w);
const double a11 =   stroke / pow(t1, 2);
const double a12 = - stroke / pow(t1, 2);
const double a21 = - stroke / pow(t2, 2);
const double a22 =   stroke / pow(t2, 2);
const double v_max = a11 * t1;
const double v_min = a21 * t2;

// Outline
vector<double> x(N);
vector<double> y(N);
vector<double> x_actual(N);
vector<double> y_actual(N);
vector<double> x_1(N);
vector<double> y_1(N);
vector<double> x_2(N);
vector<double> y_2(N);
const double theta_init = asin(e/r_base);

// Functions:
void initialize();
void calDsDtheta();
void calAlpha();
void calXRadius();
void drawSTheta();
void drawVTheta();
void drawATheta();
void drawDsDthetaS();
void drawAlphaTheta();
void drawRadiusTheta();
void drawTheoryActual();
void generateMovement();
void generateMovementS();
string formatDoubleValue(double val, int fixed);

int main() {

    initialize();
    calDsDtheta();
    calAlpha();
    calXRadius();
    drawSTheta();
    drawVTheta();
    drawATheta(); 
    drawDsDthetaS();
    drawAlphaTheta();
    drawRadiusTheta();
    drawTheoryActual(); 
    // generateMovement();
    generateMovementS(); 

    return 0;
}

string formatDoubleValue(double val, int fixed) {
    auto str = std::to_string(val);
    return str.substr(0, str.find(".") + fixed + 1);
}

void initialize() {

    for (int i = 0; i < N; i++) {
        
        theta[i] = double(360) / N * i;
        theta_pi[i] = theta[i] * PI / 180;

        // Push up
        if (theta[i] < cam_angle1) {
            if (theta[i] < cam_angle1 / 2) {
                double t = theta_pi[i] / w;
                a[i] = a11;
                v[i] = a11 * t;
                s[i] = 0.5 * a11 * pow(t, 2);
            } else {
                double t = theta_pi[i] / w - t1;
                a[i] = a12;
                v[i] = v_max + a12 * t;
                s[i] = 0.5 * stroke + v_max * t + 0.5 * a12 * pow(t, 2);
            }
        }

        // Push up stop
        else if (theta[i] >= cam_angle1 && theta[i] < cam_angle1 + stop_angle_far) {
            a[i] = 0;
            v[i] = 0;
            s[i] = stroke;
        }

        // Push down
        else if (theta[i] >= cam_angle1 + stop_angle_far && theta[i] < cam_angle1 + stop_angle_far + cam_angle2) {
            if (theta[i] < cam_angle1 + stop_angle_far + cam_angle2 / 2) {
                double t = (theta[i] - cam_angle1 - stop_angle_far) * PI / 180 / w;
                a[i] = a21;
                v[i] = a21 * t;
                s[i] = stroke + 0.5 * a21 * pow(t, 2);
            } else {
                double t = (theta[i] - cam_angle1 - stop_angle_far) * PI / 180 / w - t2;
                a[i] = a22;
                v[i] = v_min + a22 * t;
                s[i] = 0.5 * stroke + v_min * t + 0.5 * a22 * pow(t, 2);
            }
        }

        // Push down stop
        else {
            a[i] = 0;
            v[i] = 0;
            s[i] = 0;
        }

        // Check:
        // printf("theta: %.2f  a: %.3f  v: %.3f  s: %.6f\n", theta[i], a[i], v[i], s[i]);
    }
    
}

void calDsDtheta() {

    double ds = 0.0;
    double dp = 0.0;
    double ddsdp = 0.0;

    for (int i = 0; i < N; i++) {

        if (i == 0) {
            ds = s[i] - s[N-1];
            dp = theta_pi[i] - theta_pi[N-1];
        } else {
            ds = s[i] - s[i-1];
            dp = theta_pi[i] - theta_pi[i-1];
        }

        if (ds == 0.0) ds_dtheta[i] = 0.0;
        else ds_dtheta[i] = ds / dp;

        // Check:
        // printf("theta: %.2f  s: %.6f  ds_dtheta: %.6f\n", theta[i], s[i], ds_dtheta[i]);
    }

    dp = 0.0;

    for (int i = 0; i < N; i++) {

        if (i == 0) {
            ddsdp = ds_dtheta[i] - ds_dtheta[N-1];
            dp = theta_pi[i] - theta_pi[N-1];
        } else {
            ddsdp = ds_dtheta[i] - ds_dtheta[i-1];
            dp = theta_pi[i] - theta_pi[i-1];
        }

        if (ddsdp == 0.0) d2s_dtheta2[i] = 0.0;
        else d2s_dtheta2[i] = ddsdp / dp;

        // Check:
        // printf("theta: %.2f  s: %.6f  ds_dtheta: %.6f  d2s_dtheta2: %.6f\n", theta[i], s[i], ds_dtheta[i], d2s_dtheta2[i]);
    }
    
}

void calAlpha() {
    for (int i = 0; i < N; i++) {
        alpha[i] = atan2(abs(ds_dtheta[i] - e), s0 + s[i]) * 180 / PI;
        // Check:
        // printf("theta: %.3f  alpha: %.2f\n", theta[i], alpha[i]);
    }
}

void calXRadius() {
    for (int i = 0; i < N; i++) {

        x[i] = -(s0 + s[i]) * sin(theta_pi[i]) - e * cos(theta_pi[i]);
        y[i] =  (s0 + s[i]) * cos(theta_pi[i]) - e * sin(theta_pi[i]);

        // double r = sqrt(pow(x[i], 2) + pow(y[i], 2)) - r_roller;
        // x_actual[i] = x[i] * r / (r + r_roller);
        // y_actual[i] = y[i] * r / (r + r_roller);

        // double ang1 = asin(x[i]/(r+r_roller));
        // double ang2 = acos(y[i]/(r+r_roller));
        // x_actual[i] = r * sin(ang1);
        // y_actual[i] = r * cos(ang2);

        // x_actual[i] = -(s0 + s[i] - r_roller) * sin(theta_pi[i]) - e * cos(theta_pi[i]);
        // y_actual[i] =  (s0 + s[i] - r_roller) * cos(theta_pi[i]) - e * sin(theta_pi[i]);

        x_1[i] = -(s0 + s[i]) * cos(theta_pi[i]) - (ds_dtheta[i] - e) * sin(theta_pi[i]);
        y_1[i] = -(s0 + s[i]) * sin(theta_pi[i]) + (ds_dtheta[i] - e) * cos(theta_pi[i]);

        x_2[i] = -(2 * ds_dtheta[i] - e) * cos(theta_pi[i]) - (d2s_dtheta2[i] - s0 - s[i]) * sin(theta_pi[i]);
        y_2[i] =  (d2s_dtheta2[i] - s0 - s[i]) * cos(theta_pi[i]) - (2 * ds_dtheta[i] - e) * sin(theta_pi[i]);

        radius[i] = (sqrt(pow(pow(x_1[i], 2) + pow(y_1[i], 2), 3))) / (x_1[i] * y_2[i] - y_1[i] * x_2[i]);
        // radius[i] = (sqrt(pow(pow(s[i] + s0, 2) + pow(ds_dtheta[i] - e, 2), 3))) / (-(s[i] + s0) * (d2s_dtheta2[i] - s0 - s[i]) + (ds_dtheta[i] - e) * (2 * ds_dtheta[i] - e)); 

        if (abs(radius[i]) < min_radius) {
            min_radius = radius[i];
        }

        // Check:
        // printf("theta: %.3f  radius: %.3f\n", theta[i], radius[i]);
        // printf("theta: %.3f  radius: %.3f  x: %.3f  y: %.3f\n", theta[i], radius[i], x[i], y[i]);
    }

    for (int i = 0; i < N; i++) {
        
        // Previous
        double x1, y1, x2, y2;

        if (i == 0) {
            x1 = x[N - 1];
            y1 = y[N - 1];
        } else {
            x1 = x[i - 1];
            y1 = y[i - 1];
        }
        double dx1 = x[i] - x1;
        double dy1 = y[i] - y1;

        // Ahead:
        if (i == N - 1) {
            x2 = x[0];
            y2 = y[0];
        } else {
            x2 = x[i + 1];
            y2 = y[i + 1];
        }
        double dx2 = x2 - x[i];
        double dy2 = y2 - y[i];

        // Tangent vector:
        double tangent_x = (dx1 + dx2) / 2.0;
        double tangent_y = (dy1 + dy2) / 2.0;

        // Normal vector:
        double normal_x = -tangent_y;
        double normal_y =  tangent_x;

        // 计算新点的坐标
        double length = sqrt(pow(normal_x, 2) + pow(normal_y, 2));
        x_actual[i] = x[i] + r_roller * normal_x / length;
        y_actual[i] = y[i] + r_roller * normal_y / length;
        
    }

}

void drawSTheta() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 5;
    const int scale_y = 20;
    const int x_add = 300;
    const int y_add = 200;


    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(theta[i] * scale_x, s[i] * scale_y));
    
    const int x_max = int(360 * scale_x);
    const int y_max = int(stroke  * scale_y);

    // Create image and initialize it: 
    cv::Mat image_s_theta = cv::Mat::zeros(y_max + 3 * y_add, x_max + 3 * x_add, CV_8UC3);

    // Draw theta and s line:
    cv::line(image_s_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_s_theta.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_s_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_s_theta, "0", cv::Point(x_add - 50, y_max + 2 * y_add + 50), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_s_theta, "theta", cv::Point(image_s_theta.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_s_theta, "s", cv::Point(x_add - 15, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_s_theta, cv::Point(x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // Up stop start:
    cv::line(image_s_theta, cv::Point(x_add + cam_angle1 * scale_x, 2 * y_add), cv::Point(x_add + cam_angle1 * scale_x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_s_theta, cv::Point(x_add + cam_angle1 * scale_x, 2 * y_add), cv::Point(x_add, 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_s_theta, formatDoubleValue(cam_angle1, 1), cv::Point(x_add + cam_angle1 * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_s_theta, formatDoubleValue(stroke, 1), cv::Point(x_add - 140, 2 * y_add + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // Up stop finished Down start:
    cv::line(image_s_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x, 2 * y_add), cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_s_theta, formatDoubleValue(cam_angle1 + stop_angle_far, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // Down stop start:
    cv::putText(image_s_theta, formatDoubleValue(cam_angle1 + stop_angle_far + cam_angle2, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2) * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    
    // 360
    cv::putText(image_s_theta, "360.0", cv::Point(x_add + x_max - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    
    // Show images:
    // cv::imshow("image_s_theta", image_s_theta);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_s_theta.tiff", image_s_theta);
    cv::imwrite("../images/JPEG/image_s_theta.jpeg", image_s_theta);

}

void drawVTheta() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 5;
    const int scale_y = 2;
    const int x_add = 300;
    const int y_add = 200;


    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(theta[i] * scale_x, v[i] * scale_y));

    const int x_max = int(360 * scale_x);
    const int y_max = int(v_max * scale_y);
    const int y_min = int(v_min * scale_y);

    // Create image and initialize it: 
    cv::Mat image_v_theta = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max + 3 * x_add, CV_8UC3);

    // Draw theta and v line:
    cv::line(image_v_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_v_theta.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_v_theta, cv::Point(x_add, y_max - y_min + 3 * y_add), cv::Point(x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_v_theta, "0", cv::Point(x_add - 50, y_max + 2 * y_add + 50), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_v_theta, "theta", cv::Point(image_v_theta.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_v_theta, "v", cv::Point(x_add - 15, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_v_theta, cv::Point(x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // v_max
    cv::circle(image_v_theta, cv::Point(x_add + cam_angle1 / 2 * scale_x, 2 * y_add), 2, cv::Scalar(10,215,255), -1);
    cv::line(image_v_theta, cv::Point(x_add + cam_angle1 / 2 * scale_x, 2 * y_add), cv::Point(x_add + cam_angle1 / 2 * scale_x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_v_theta, cv::Point(x_add + cam_angle1 / 2 * scale_x, 2 * y_add), cv::Point(x_add, 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_v_theta, formatDoubleValue(cam_angle1 / 2, 1), cv::Point(x_add + cam_angle1 / 2 * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_v_theta, formatDoubleValue(v_max, 1), cv::Point(x_add - 170, 2 * y_add + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    cv::putText(image_v_theta, formatDoubleValue(cam_angle1, 1), cv::Point(x_add + cam_angle1 * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // v_min
    cv::circle(image_v_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x,  y_max - y_min + 2 * y_add), 2, cv::Scalar(10,215,255), -1);
    cv::line(image_v_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x,  y_max - y_min + 2 * y_add), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x,  y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_v_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x,  y_max - y_min + 2 * y_add), cv::Point(x_add,  y_max - y_min + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_v_theta, formatDoubleValue(cam_angle1 + stop_angle_far + cam_angle2 / 2, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x - 60, y_max + 2 * y_add - 27), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_v_theta, formatDoubleValue(v_min, 1), cv::Point(x_add - 200, y_max - y_min + 2 * y_add + 20), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    cv::putText(image_v_theta, formatDoubleValue(cam_angle1 + stop_angle_far + cam_angle2, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2) * scale_x - 60, y_max + 2 * y_add - 27), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);


    // 360
    cv::putText(image_v_theta, "360.0", cv::Point(x_add + x_max - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    
    // Show images:
    // cv::imshow("image_v_theta", image_v_theta);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_v_theta.tiff", image_v_theta);
    cv::imwrite("../images/JPEG/image_v_theta.jpeg", image_v_theta);

}

void drawATheta() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 5;
    const double scale_y = 0.2;
    const int x_add = 300;
    const int y_add = 200;

    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(theta[i] * scale_x, a[i] * scale_y));

    const int x_max = int(360 * scale_x);
    const int y_max = int(max(a11, a22) * scale_y);
    const int y_min = int(min(a12, a21) * scale_y);

    // Create image and initialize it: 
    cv::Mat image_a_theta = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max + 3 * x_add, CV_8UC3);

    // Draw theta and v line:
    cv::line(image_a_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_a_theta.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_a_theta, cv::Point(x_add, y_max - y_min + 3 * y_add), cv::Point(x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_a_theta, "0", cv::Point(x_add - 50, y_max + 2 * y_add + 50), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_a_theta, "theta", cv::Point(image_a_theta.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_a_theta, "a", cv::Point(x_add - 20, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_a_theta, cv::Point(x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // a11
    cv::line(image_a_theta, cv::Point(x_add + cam_angle1 / 2 * scale_x, y_max + 2 * y_add - a11 * scale_y), cv::Point(x_add + cam_angle1 / 2 * scale_x, y_max + 2 * y_add - a12 * scale_y), cv::Scalar(230,216,173), 1);
    cv::putText(image_a_theta, formatDoubleValue(cam_angle1 / 2, 1), cv::Point(x_add + cam_angle1 / 2 * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_a_theta, formatDoubleValue(a11, 1), cv::Point(x_add - 200, y_max + 2 * y_add - a11 * scale_y + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // a12
    cv::line(image_a_theta, cv::Point(x_add + cam_angle1 / 2 * scale_x, y_max + 2 * y_add - a12 * scale_y), cv::Point(x_add, y_max + 2 * y_add - a12 *scale_y), cv::Scalar(230,216,173), 1);
    cv::line(image_a_theta, cv::Point(x_add + cam_angle1 * scale_x, y_max + 2 * y_add - a12 * scale_y), cv::Point(x_add + cam_angle1 * scale_x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_a_theta, formatDoubleValue(cam_angle1, 1), cv::Point(x_add + cam_angle1 * scale_x - 60, y_max + 2 * y_add - 27), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_a_theta, formatDoubleValue(a12, 1), cv::Point(x_add - 225, y_max + 2 * y_add - a12 * scale_y + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // a21
    cv::line(image_a_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x,  y_max - a21 * scale_y + 2 * y_add), cv::Point(x_add, y_max - a21 * scale_y + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_a_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x,  y_max - a21 * scale_y + 2 * y_add), cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x,  y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_a_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x,  y_max - a21 * scale_y + 2 * y_add), cv::Point(x_add + (cam_angle1 + stop_angle_far+ cam_angle2 / 2) * scale_x,  y_max + 2 * y_add - a22 * scale_y), cv::Scalar(230,216,173), 1);
    cv::putText(image_a_theta, formatDoubleValue(cam_angle1 + stop_angle_far, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far) * scale_x - 60, y_max + 2 * y_add - 27), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_a_theta, formatDoubleValue(cam_angle1 + stop_angle_far + cam_angle2 / 2, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x - 60, y_max + 2 * y_add - 27), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_a_theta, formatDoubleValue(a21, 1), cv::Point(x_add - 235, y_max - a21 * scale_y + 2 * y_add + 20), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // a22
    cv::line(image_a_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2 / 2) * scale_x, y_max + 2 * y_add - a22 * scale_y), cv::Point(x_add, y_max + 2 * y_add - a22 *scale_y), cv::Scalar(230,216,173), 1);
    cv::line(image_a_theta, cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2) * scale_x, y_max + 2 * y_add - a22 * scale_y), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2) * scale_x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_a_theta, formatDoubleValue(cam_angle1 + stop_angle_far + cam_angle2, 1), cv::Point(x_add + (cam_angle1 + stop_angle_far + cam_angle2) * scale_x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_a_theta, formatDoubleValue(a22, 1), cv::Point(x_add - 200, y_max + 2 * y_add - a22 * scale_y + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // 360
    cv::putText(image_a_theta, "360.0", cv::Point(x_add + x_max - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    
    // Show images:
    // cv::imshow("image_a_theta", image_a_theta);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_a_theta.tiff", image_a_theta);
    cv::imwrite("../images/JPEG/image_a_theta.jpeg", image_a_theta);

}

void drawDsDthetaS() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 10;
    const int scale_y = 10;
    const int x_add = 400;
    const int y_add = 200;

    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(ds_dtheta[i] * scale_x, s[i] * scale_y));

    double x_max_temp = 0.0;
    int x_max_index = 0;
    double x_min_temp = 0.0;
    int x_min_index = 0;
    double x_temp = 0.1;
    int x_index = 0;

    for (int i = 0; i < N; i++) {

        if (theta[i] > cam_angle1 + stop_angle_far +  cam_angle2 / 2 && theta[i] <= cam_angle1 + stop_angle_far + cam_angle2) {
            double x1, y1, x2, y2;
            
            x1 = ds_dtheta[i - 1];
            y1 = s[i - 1];
            
            double dx1 = ds_dtheta[i] - x1;
            double dy1 = s[i] - y1;

            if (i == N - 1) {
                x2 = ds_dtheta[0];
                y2 = s[0];
            } else {
                x2 = ds_dtheta[i + 1];
                y2 = s[i + 1];
            }
            double dx2 = x2 - ds_dtheta[i];
            double dy2 = y2 - s[i];

            // Tangent vector:
            double tangent_x = dx1 + dx2;
            double tangent_y = dy1 + dy2;
            double temp = abs(atan2(tangent_x, -tangent_y) - pressure_angle2 * PI / 180);

            printf("%.3f:  %.3f\n", theta[i], temp);

            if (temp < x_temp) {
                x_temp = temp;
                x_index = i;
            }
        }
        
        if (ds_dtheta[i] > x_max_temp) {
            x_max_temp = ds_dtheta[i];
            x_max_index = i;
        }
        if (ds_dtheta[i] < x_min_temp) {
            x_min_temp = ds_dtheta[i];
            x_min_index = i;
        }
    }

    printf("\n%.3f:  x_index: %d\n", theta[x_index], x_index);

    const int x_max = int(x_max_temp * scale_x);
    const int x_min = int(x_min_temp * scale_x);
    const int y_max = int(stroke * scale_y);
    const int y_min = int(-s0 * scale_y);

    // Create image and initialize it: 
    cv::Mat image_ds_dtheta_s = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max - x_min + 4 * x_add, CV_8UC3);

    // Draw theta and v line:
    cv::line(image_ds_dtheta_s, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_ds_dtheta_s.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add, y_max - y_min + 3 * y_add), cv::Point(-x_min + 2 * x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_ds_dtheta_s, "0", cv::Point(-x_min + 2 * x_add - 50, y_max + 2 * y_add + 70), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_ds_dtheta_s, "ds/dtheta", cv::Point(image_ds_dtheta_s.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_ds_dtheta_s, "s", cv::Point(-x_min + 2 * x_add - 15, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // x_max
    cv::circle(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), 2, cv::Scalar(10,215,255), -1);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Scalar(230,216,173), 1);
    cv::putText(image_ds_dtheta_s, formatDoubleValue(x_max_temp, 2), cv::Point(-x_min + 2 * x_add + position[x_max_index].x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_ds_dtheta_s, formatDoubleValue(stroke / 2, 1), cv::Point(-x_min + 2 * x_add + 30, y_max + 2 * y_add - stroke / 2 * scale_y - 30), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // x_min
    cv::circle(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), 2, cv::Scalar(10,215,255), -1);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Scalar(230,216,173), 1);
    cv::putText(image_ds_dtheta_s, formatDoubleValue(x_min_temp, 2), cv::Point(-x_min + 2 * x_add + position[x_min_index].x - 100, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    
    // s0 e pressure_angle1 and pressure_angle2:
    cv::circle(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + e * scale_x, y_max + 2 * y_add + s0 * scale_y), 8, cv::Scalar(10,215,255), -1);
    cv::putText(image_ds_dtheta_s, "(e, -s0)", cv::Point(-x_min + 2 * x_add + e * scale_x + 15, y_max + 2 * y_add + s0 * scale_y + 70), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    const double k1 =  1 / tan(pressure_angle1 * PI / 180);
    const double k2 = -1 / tan(pressure_angle1 * PI / 180);
    const double k3 =  1 / tan(pressure_angle2 * PI / 180);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), cv::Point(-x_min + 2 * x_add + y_min / k2, y_max - y_min + 2 * y_add), cv::Scalar(203,192,255), 2);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - stroke / 2 * scale_y), cv::Point(-x_min + 2 * x_add + position[x_max_index].x - (-y_min + stroke / 2 * scale_y) / k1, y_max - y_min + 2 * y_add), cv::Scalar(203,192,255), 2);
    cv::line(image_ds_dtheta_s, cv::Point(-x_min + 2 * x_add + position[x_index].x, y_max + 2 * y_add - position[x_index].y), cv::Point(-x_min + 2 * x_add + position[x_min_index].x + (-y_min + position[x_index].y) / k3, y_max - y_min + 2 * y_add), cv::Scalar(203,192,255), 2);

    // Show images:
    // cv::imshow("image_ds_dtheta_s", image_ds_dtheta_s);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_ds_dtheta_s.tiff", image_ds_dtheta_s);
    cv::imwrite("../images/JPEG/image_ds_dtheta_s.jpeg", image_ds_dtheta_s);

}

void drawAlphaTheta() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 5;
    const int scale_y = 20;
    const int x_add = 300;
    const int y_add = 200;

    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(theta[i] * scale_x, alpha[i] * scale_y));

    double y_up = 0.0;
    int y_up_index = 0;
    double y_down = 0.0;
    int y_down_index = 0;

    for (int i = 0; i < N; i++) {
        if (theta[i] < cam_angle1) {
            if (alpha[i] > y_up) {
                y_up = alpha[i];
                y_up_index = i;
            }
        } else {
            if (alpha[i] > y_down) {
                y_down = alpha[i];
                y_down_index = i;
            }
        }
    }

    const int x_max = int(360 * scale_x);
    const int y_max = int(max(y_up, y_down) * scale_y);
    

    // Create image and initialize it: 
    cv::Mat image_alpha_theta = cv::Mat::zeros(y_max + 3 * y_add, x_max + 3 * x_add, CV_8UC3);

    // Draw theta and alpha line:
    cv::line(image_alpha_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_alpha_theta.cols - x_add, y_max +2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_alpha_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_alpha_theta, "0", cv::Point(x_add - 50, y_max + 2 * y_add + 50), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_alpha_theta, "theta", cv::Point(image_alpha_theta.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_alpha_theta, "alpha", cv::Point(x_add - 70, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_alpha_theta, cv::Point(x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // alpha_up
    // cv::line(image_alpha_theta, cv::Point(x_add + int(position[y_up_index].x), y_max + 2 * y_add - int(position[y_up_index].y)), cv::Point(x_add + int(position[y_up_index].x), y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    // cv::line(image_alpha_theta, cv::Point(x_add + int(position[y_up_index].x), y_max + 2 * y_add - int(position[y_up_index].y)), cv::Point(x_add, y_max + 2 * y_add - int(position[y_up_index].y)), cv::Scalar(230,216,173), 1);
    // cv::putText(image_alpha_theta, formatDoubleValue(theta[y_up_index], 1), cv::Point(x_add + int(position[y_up_index].x) - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_alpha_theta, formatDoubleValue(y_up, 1), cv::Point(x_add - 150, y_max + 2 * y_add - int(position[y_up_index].y) + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // alpha_down
    cv::line(image_alpha_theta, cv::Point(x_add + int(position[y_down_index].x), y_max + 2 * y_add - int(position[y_down_index].y)), cv::Point(x_add + int(position[y_down_index].x), y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_alpha_theta, cv::Point(x_add + int(position[y_down_index].x), y_max + 2 * y_add - int(position[y_down_index].y)), cv::Point(x_add, y_max + 2 * y_add - int(position[y_down_index].y)), cv::Scalar(230,216,173), 1);
    cv::putText(image_alpha_theta, formatDoubleValue(theta[y_down_index], 1), cv::Point(x_add + int(position[y_down_index].x) - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_alpha_theta, formatDoubleValue(y_down, 1), cv::Point(x_add - 150, y_max + 2 * y_add - int(position[y_down_index].y) + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // 360
    cv::putText(image_alpha_theta, "360.0", cv::Point(x_add + x_max - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    
    // Show images:
    // cv::imshow("image_alpha_theta", image_alpha_theta);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_alpha_theta.tiff", image_alpha_theta);
    cv::imwrite("../images/JPEG/image_alpha_theta.jpeg", image_alpha_theta);

}

void drawRadiusTheta() {

    vector<Point2f> position;

    // Procession of x and y:
    const int scale_x = 5;
    const int scale_y = 2;
    const int x_add = 300;
    const int y_add = 200;

    double y_max_temp = 0.0;
    int y_max_index = 0;
    double y_min_temp = 100.0;
    int y_min_index = 0;

    for (int i = 0; i < N; i++) {
        if (radius[i] > y_max_temp) {
            y_max_temp = radius[i];
            y_max_index = i;
        }
        if (radius[i] < y_min_temp) {
            y_min_temp = radius[i];
            y_min_index = i;
        }
    }

    const int x_max = int(360 * scale_x);
    const int y_max = int(y_max_temp * scale_y);
    const int y_min = int(y_min_temp * scale_y);


    // Put all points into vector:
    for (int i = 0; i < N; i++) 
        position.push_back(Point2f(theta[i] * scale_x, radius[i] * scale_y));


    // Create image and initialize it: 
    cv::Mat image_radius_theta = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max + 3 * x_add, CV_8UC3);

    // Draw theta and v line:
    cv::line(image_radius_theta, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_radius_theta.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 1);
    cv::line(image_radius_theta, cv::Point(x_add, y_max - y_min + 3 * y_add), cv::Point(x_add, y_add), cv::Scalar(255,255,255), 1);
    
    // Draw initial point(0) and texts:
    cv::putText(image_radius_theta, "0", cv::Point(x_add - 50, y_max + 2 * y_add + 50), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(255,255,255), 1, 6);
    cv::putText(image_radius_theta, "theta", cv::Point(image_radius_theta.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    cv::putText(image_radius_theta, "radius", cv::Point(x_add - 70, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(144,238,144), 1, 6);
    
    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::circle(image_radius_theta, cv::Point(x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 2, cv::Scalar(10,215,255), -1);
    }

    // radius_max
    cv::line(image_radius_theta, cv::Point(x_add + position[y_max_index].x, 2 * y_add), cv::Point(x_add + position[y_max_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_radius_theta, cv::Point(x_add + position[y_max_index].x, 2 * y_add), cv::Point(x_add, 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::putText(image_radius_theta, formatDoubleValue(theta[y_max_index], 1), cv::Point(x_add + position[y_max_index].x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_radius_theta, formatDoubleValue(y_max_temp, 1), cv::Point(x_add - 170, 2 * y_add + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // radius_min
    cv::line(image_radius_theta, cv::Point(x_add + position[y_min_index].x, y_max + 2 * y_add - position[y_min_index].y), cv::Point(x_add + position[y_min_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    cv::line(image_radius_theta, cv::Point(x_add + position[y_min_index].x, y_max + 2 * y_add - position[y_min_index].y), cv::Point(x_add, y_max + 2 * y_add - position[y_min_index].y), cv::Scalar(230,216,173), 1);
    cv::putText(image_radius_theta, formatDoubleValue(theta[y_min_index], 1), cv::Point(x_add + position[y_min_index].x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_radius_theta, formatDoubleValue(y_min_temp, 1), cv::Point(x_add - 150, y_max + 2 * y_add - position[y_min_index].y + 15), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    // 360
    cv::putText(image_radius_theta, "360.0", cv::Point(x_add + x_max - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);

    
    // Show images:
    // cv::imshow("image_radius_theta", image_radius_theta);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_radius_theta.tiff", image_radius_theta);
    cv::imwrite("../images/JPEG/image_radius_theta.jpeg", image_radius_theta);

}

void drawTheoryActual() {

    vector<Point2f> position;
    vector<Point2f> position_actual;

    // Procession of x and y:
    const int scale_x = 8;
    const int scale_y = 8;
    const int x_add = 300;
    const int y_add = 300;

    // Put all points into vector:
    for (int i = 0; i < N; i++) {
        position.push_back(Point2f(x[i] * scale_x, y[i] * scale_y));
        position_actual.push_back(Point2f(x_actual[i] * scale_x, y_actual[i] * scale_y));
    }

    double x_max_temp = 0.0;
    int x_max_index = 0;
    double x_min_temp = 0.0;
    int x_min_index = 0;

    double y_max_temp = 0.0;
    int y_max_index = 0;
    double y_min_temp = 0.0;
    int y_min_index = 0;

    for (int i = 0; i < N; i++) {
        if (x[i] > x_max_temp) {
            x_max_temp = x[i];
            x_max_index = i;
        }
        if (x[i] < x_min_temp) {
            x_min_temp = x[i];
            x_min_index = i;
        }
        if (y[i] > y_max_temp) {
            y_max_temp = y[i];
            y_max_index = i;
        }
        if (y[i] < y_min_temp) {
            y_min_temp = y[i];
            y_min_index = i;
        }
    }

    const int x_max = int(x_max_temp * scale_x);
    const int x_min = int(x_min_temp * scale_x);
    const int y_max = int(y_max_temp * scale_y);
    const int y_min = int(y_min_temp * scale_y);

    // Create image and initialize it: 
    cv::Mat image_theory_actual = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max - x_min + 4 * x_add, CV_8UC3);

    // Draw theta and v line:
    cv::line(image_theory_actual, cv::Point(x_add, y_max + 2 * y_add), cv::Point(image_theory_actual.cols - x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 3);
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add, y_max - y_min + 3 * y_add), cv::Point(-x_min + 2 * x_add, y_add), cv::Scalar(255,255,255), 3);
    
    // Draw initial point(0) and texts:
    cv::putText(image_theory_actual, "0", cv::Point(-x_min + 2 * x_add - 70, y_max + 2 * y_add + 85), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(255,255,255), 2, 6);
    cv::putText(image_theory_actual, "x", cv::Point(image_theory_actual.cols - x_add + 40, y_max + 2 * y_add + 18), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
    cv::putText(image_theory_actual, "y", cv::Point(-x_min + 2 * x_add - 20, y_add - 60), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
    
    // R_base and r_e
    cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), r_base * scale_x, cv::Scalar(144,238,144), 3);
    cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), e * scale_x, cv::Scalar(144,238,144), 3);
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add), cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y[0] * scale_y), cv::Scalar(230,216,173), 1);
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add - e * scale_x, y_add), cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y[0] * scale_y), cv::Scalar(20,80,255), 10);
    cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y[0] * scale_y), 2.5 * scale_x, cv::Scalar(20,80,255), -1);
    cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y[0] * scale_y), r_roller * scale_x, cv::Scalar(20,80,255), 12);

    // Draw points on mat:
    for (int i = 0; i < position.size(); i++) {
        cv::Point2f pos = position[i];
        cv::Point2f pos_actual = position_actual[i];
        cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add + int(pos.x), y_max + 2 * y_add - int(pos.y)), 4, cv::Scalar(10,215,255), -1);
        cv::circle(image_theory_actual, cv::Point(-x_min + 2 * x_add + int(pos_actual.x), y_max + 2 * y_add - int(pos_actual.y)), 4, cv::Scalar(147,20,255), -1);
    }

    // x_max
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - position[x_max_index].y), cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 2);
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[x_max_index].x, y_max + 2 * y_add - position[x_max_index].y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - position[x_max_index].y), cv::Scalar(230,216,173), 2);
    cv::putText(image_theory_actual, formatDoubleValue(x_max_temp, 1), cv::Point(-x_min + 2 * x_add + position[x_max_index].x - 30, y_max + 2 * y_add - 50), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
    cv::putText(image_theory_actual, formatDoubleValue(position[x_max_index].y / scale_y, 1), cv::Point(-x_min + 2 * x_add + 30, y_max + 2 * y_add - position[x_max_index].y + 80), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);

    // x_min
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add - position[x_min_index].y), cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 2);
    cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[x_min_index].x, y_max + 2 * y_add - position[x_min_index].y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - position[x_min_index].y), cv::Scalar(230,216,173), 2);
    cv::putText(image_theory_actual, formatDoubleValue(x_min_temp, 1), cv::Point(-x_min + 2 * x_add + position[x_min_index].x - 200, y_max + 2 * y_add - 60), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
    cv::putText(image_theory_actual, formatDoubleValue(position[x_min_index].y / scale_y, 1), cv::Point(-x_min + 2 * x_add - 220, y_max + 2 * y_add - position[x_min_index].y + 80), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
    
    // y_max
    // cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[y_max_index].x, y_max + 2 * y_add - position[y_max_index].y), cv::Point(-x_min + 2 * x_add + position[y_max_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    // cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[y_max_index].x, y_max + 2 * y_add - position[y_max_index].y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - position[y_max_index].y), cv::Scalar(230,216,173), 1);
    // cv::putText(image_theory_actual, formatDoubleValue(position[x_max_index].x / scale_x, 1), cv::Point(-x_min + 2 * x_add + position[y_max_index].x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_theory_actual, formatDoubleValue(y_max_temp, 1), cv::Point(-x_min + 2 * x_add + 35, y_max + 2 * y_add - position[y_max_index].y - 35), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);

    // y_min
    // cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[y_min_index].x, y_max + 2 * y_add - position[y_min_index].y), cv::Point(-x_min + 2 * x_add + position[y_min_index].x, y_max + 2 * y_add), cv::Scalar(230,216,173), 1);
    // cv::line(image_theory_actual, cv::Point(-x_min + 2 * x_add + position[y_min_index].x, y_max + 2 * y_add - position[y_min_index].y), cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add - position[y_min_index].y), cv::Scalar(230,216,173), 1);
    // cv::putText(image_theory_actual, formatDoubleValue(position[y_min_index].x / scale_x, 1), cv::Point(-x_min + 2 * x_add + position[y_min_index].x - 60, y_max + 2 * y_add + 60), cv::FONT_HERSHEY_SIMPLEX, 1.5, cv::Scalar(0,69,255), 1, 6);
    cv::putText(image_theory_actual, formatDoubleValue(y_min_temp, 1), cv::Point(-x_min + 2 * x_add - 260, y_max + 2 * y_add - position[y_min_index].y + 90), cv::FONT_HERSHEY_SIMPLEX, 2, cv::Scalar(0,69,255), 2, 6);
    
    

    // Show images:
    // cv::imshow("image_theory_actual", image_theory_actual);
    // cv::waitKey();

    // Save image:
    cv::imwrite("../images/TIFF/image_theory_actual.tiff", image_theory_actual);
    cv::imwrite("../images/JPEG/image_theory_actual.jpeg", image_theory_actual);

}

void generateMovement() {

    // Initialize MP4V codec
    int fourcc = cv::VideoWriter::fourcc('m', 'p', '4', 'v'); 
    double fps = 60.0;

    // Create images:
    vector<Mat> images;

    // Procession of x and y:
    const int scale_x = 8;
    const int scale_y = 8;
    const int x_add = 500;
    const int y_add = 500;

    double x_max_temp = 0.0;
    int x_max_index = 0;
    double x_min_temp = 0.0;
    int x_min_index = 0;

    double y_max_temp = 0.0;
    int y_max_index = 0;
    double y_min_temp = 0.0;
    int y_min_index = 0;

    for (int i = 0; i < N; i++) {
        if (x[i] > x_max_temp) {
            x_max_temp = x[i];
            x_max_index = i;
        }
        if (x[i] < x_min_temp) {
            x_min_temp = x[i];
            x_min_index = i;
        }
        if (y[i] > y_max_temp) {
            y_max_temp = y[i];
            y_max_index = i;
        }
        if (y[i] < y_min_temp) {
            y_min_temp = y[i];
            y_min_index = i;
        }
    }

    const int x_max = int(x_max_temp * scale_x);
    const int x_min = int(x_min_temp * scale_x);
    const int y_max = int(y_max_temp * scale_y);
    const int y_min = int(y_min_temp * scale_y);

    // Create Vedio Writer:
    cv::Mat image_temp = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max - x_min + 4 * x_add, CV_8UC3);
    cv::VideoWriter video("../images/VIDEO/Movements.mp4", fourcc, fps, image_temp.size());
    
    // Generate video:
    for (int i = 0; i < N; i++) {

        double rotation = theta_pi[i];

        // Create image:
        cv::Mat image = cv::Mat::zeros(int(y_max-y_min)+4*y_add, int(x_max-x_min)+4*x_add, CV_8UC3);

        // Put name and author
        cv::putText(image, "Movements", cv::Point(image.cols-2*x_add, y_add-20), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(230,216,173), 3, 6);
        cv::putText(image, "Author: ShiJiaxiao", cv::Point(image.cols-2*x_add-120, y_add+120), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(230,216,173), 3, 6);

        {   // Draw initial point(0) and texts:
            cv::line(image, cv::Point(x_add, int(y_max)+2*y_add), cv::Point(image.cols-x_add, int(y_max)+2*y_add), cv::Scalar(255,255,255), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add, y_add), cv::Point(int(-x_min)+2*x_add, image.rows-y_add), cv::Scalar(255,255,255), 2);
            cv::putText(image, "0", cv::Point(int(-x_min)+2*x_add-80, int(y_max)+2*y_add+80), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(255,255,255), 2, 6);
            cv::putText(image, "x", cv::Point(image.cols-x_add+40, int(y_max)+2*y_add+20), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
            cv::putText(image, "y", cv::Point(int(-x_min)+2*x_add-20, y_add-60), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
        }
        {   // Draw r_base r_e and roller:
            cv::circle(image, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), r_base * scale_x, cv::Scalar(144,238,144), 3);
            cv::circle(image, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), e * scale_x, cv::Scalar(144,238,144), 3);
            double x_temp =  x[i] * cos(rotation) + y[i] * sin(rotation);
            double y_temp = -x[i] * sin(rotation) + y[i] * cos(rotation);
            cv::line(image, cv::Point(int(-x_min)+2*x_add - e * scale_x, y_max + 2 * y_add), cv::Point(int(-x_min)+2*x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), cv::Scalar(230,216,173), 2);
            cv::line(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_add), cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), cv::Scalar(20,80,255), 10);
            cv::circle(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), 2.5 * scale_x, cv::Scalar(20,80,255), -1);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_temp * scale_x, y_max + 2 * y_add - y_temp *scale_y), r_roller * scale_x, cv::Scalar(20,80,255), 12);
        }

        // Draw points
        double x_r = 0.0;
        double y_r = 0.0;
        double x_ar = 0.0;
        double y_ar = 0.0;
        
        for (int j = 0; j < N; j++) {
            x_r =  x[j] * cos(rotation) + y[j] * sin(rotation);
            y_r = -x[j] * sin(rotation) + y[j] * cos(rotation);
            x_ar =  x_actual[j] * cos(rotation) + y_actual[j] * sin(rotation);
            y_ar = -x_actual[j] * sin(rotation) + y_actual[j] * cos(rotation);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_r * scale_x, y_max + 2 * y_add - y_r *scale_y), 4, cv::Scalar(10,215,255), -1);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_ar * scale_x, y_max + 2 * y_add - y_ar *scale_y), 4, cv::Scalar(147,20,255), -1);
        }

        // Add image to video:
        images.push_back(image);
        video.write(image);

    }

    for (int i = 0 ; i < 300; i++) {
        video.write(images[images.size()-1]);
    }

    // Release VideoWriter:
    video.release();
    

}

void generateMovementS() {

    // Initialize MP4V codec
    int fourcc = cv::VideoWriter::fourcc('m', 'p', '4', 'v'); 
    double fps = 60.0;

    // Create images:
    vector<Mat> images;

    // Procession of x and y:
    const int scale_x = 5;
    const int scale_y = 5;
    const int x_add = 500;
    const int y_add = 500;

    double x_max_temp = 0.0;
    int x_max_index = 0;
    double x_min_temp = 0.0;
    int x_min_index = 0;

    double y_max_temp = 0.0;
    int y_max_index = 0;
    double y_min_temp = 0.0;
    int y_min_index = 0;

    for (int i = 0; i < N; i++) {
        if (x[i] > x_max_temp) {
            x_max_temp = x[i];
            x_max_index = i;
        }
        if (x[i] < x_min_temp) {
            x_min_temp = x[i];
            x_min_index = i;
        }
        if (y[i] > y_max_temp) {
            y_max_temp = y[i];
            y_max_index = i;
        }
        if (y[i] < y_min_temp) {
            y_min_temp = y[i];
            y_min_index = i;
        }
    }

    const int x_max = int(x_max_temp * scale_x);
    const int x_min = int(x_min_temp * scale_x);
    const int y_max = int(y_max_temp * scale_y);
    const int y_min = int(y_min_temp * scale_y);

    // Create Vedio Writer:
    cv::Mat image_temp = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max - x_min + 7 * x_add + 360 * scale_x, CV_8UC3);
    cv::VideoWriter video("../images/VIDEO/Movements_S.mp4", fourcc, fps, image_temp.size());
    
    // Generate video:
    for (int i = 0; i < N; i++) {

        double rotation = theta_pi[i];
        double x_temp =  x[i] * cos(rotation) + y[i] * sin(rotation);
        double y_temp = -x[i] * sin(rotation) + y[i] * cos(rotation);

        // Create image:
        cv::Mat image = cv::Mat::zeros(y_max - y_min + 4 * y_add, x_max - x_min + 7 * x_add + 360 * scale_x, CV_8UC3);

        // Put name and author
        cv::putText(image, "Movements_S", cv::Point(x_max - x_min + 5 * x_add + 500, y_max + 2 * y_add + 300), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(230,216,173), 3, 6);
        cv::putText(image, "Author: ShiJiaxiao", cv::Point(x_max - x_min + 5 * x_add + 400, y_max + 2 * y_add + 500), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(230,216,173), 3, 6);

        {   // Draw initial point(0) and texts:
            cv::line(image, cv::Point(x_add, y_max + 2 * y_add), cv::Point(x_max - x_min + 3 * x_add, y_max + 2 * y_add), cv::Scalar(255,255,255), 2);
            cv::line(image, cv::Point(int(-x_min)+2*x_add, y_add), cv::Point(int(-x_min)+2*x_add, image.rows-y_add), cv::Scalar(255,255,255), 2);
            cv::putText(image, "0", cv::Point(int(-x_min)+2*x_add-70, int(y_max)+2*y_add+70), cv::FONT_HERSHEY_SIMPLEX, 2.2, cv::Scalar(255,255,255), 2, 6);
            cv::putText(image, "x", cv::Point(x_max - x_min + 3 * x_add + 40, int(y_max)+2*y_add+20), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
            cv::putText(image, "y", cv::Point(int(-x_min)+2*x_add-20, y_add-60), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
        }
        {   // Draw r_base r_e and roller:
            cv::circle(image, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), r_base * scale_x, cv::Scalar(144,238,144), 3);
            cv::circle(image, cv::Point(-x_min + 2 * x_add, y_max + 2 * y_add), e * scale_x, cv::Scalar(144,238,144), 3);
            cv::line(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add), cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), cv::Scalar(230,216,173), 2);
            cv::line(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_add), cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), cv::Scalar(20,80,255), 10);
            cv::circle(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), 2.5 * scale_x, cv::Scalar(20,80,255), -1);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_temp * scale_x, y_max + 2 * y_add - y_temp *scale_y), r_roller * scale_x, cv::Scalar(20,80,255), 12);
        }

        // Draw points
        double x_r = 0.0;
        double y_r = 0.0;
        double x_ar = 0.0;
        double y_ar = 0.0;
        
        for (int j = 0; j < N; j++) {
            x_r =  x[j] * cos(rotation) + y[j] * sin(rotation);
            y_r = -x[j] * sin(rotation) + y[j] * cos(rotation);
            x_ar =  x_actual[j] * cos(rotation) + y_actual[j] * sin(rotation);
            y_ar = -x_actual[j] * sin(rotation) + y_actual[j] * cos(rotation);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_r * scale_x, y_max + 2 * y_add - y_r * scale_y), 4, cv::Scalar(10,215,255), -1);
            cv::circle(image, cv::Point(-x_min + 2 * x_add + x_ar * scale_x, y_max + 2 * y_add - y_ar * scale_y), 4, cv::Scalar(147,20,255), -1);
        }


        // Draw S:
        cv::line(image, cv::Point(x_max - x_min + 5 * x_add, y_max + 2 * y_add - s0 * scale_y), cv::Point(image.cols - x_add, y_max + 2 * y_add - s0 * scale_y), cv::Scalar(255,255,255), 2);
        cv::line(image, cv::Point(x_max - x_min + 5 * x_add, y_max + 2 * y_add - s0 * scale_y), cv::Point(x_max - x_min + 5 * x_add, y_add), cv::Scalar(255,255,255), 2);
        
        // Draw initial point(0) and texts:
        cv::putText(image, "0", cv::Point(x_max - x_min + 5 * x_add - 70, y_max + 2 * y_add - s0 * scale_y + 70), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(255,255,255), 2, 6);
        cv::putText(image, "theta", cv::Point(image.cols - x_add + 40, y_max + 2 * y_add - s0 * scale_y + 22), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
        cv::putText(image, "s", cv::Point(x_max - x_min + 5 * x_add - 25, y_add - 40), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(144,238,144), 2, 6);
        
        // Draw points on mat:
        for (int j = 0; j < i; j++) {
            cv::circle(image, cv::Point(x_max - x_min + 5 * x_add + theta[j] * scale_x, y_max + 2 * y_add - (s0+s[j]) * scale_y), 4, cv::Scalar(203,192,255), -1);
        }

        // S point
        cv::line(image, cv::Point(-x_min + 2 * x_add - e * scale_x, y_max + 2 * y_add - y_temp * scale_y), cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x, y_max + 2 * y_add - (s0+s[i]) * scale_y), cv::Scalar(230,216,173), 2);
        cv::line(image, cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x, y_max + 2 * y_add - (s0+s[i]) * scale_y), cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x, y_max + 2 * y_add - s0 * scale_y), cv::Scalar(230,216,173), 2);
        if (!strcmp(formatDoubleValue(theta[i], 1).c_str(), "0.3"))
            cv::putText(image, "0.0", cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x + 25, y_max + 2 * y_add - s0 * scale_y + 125), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(0,69,255), 2, 6);
        else if (!strcmp(formatDoubleValue(theta[i], 1).c_str(), "359.6"))
            cv::putText(image, "360.0", cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x + 25, y_max + 2 * y_add - s0 * scale_y + 125), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(0,69,255), 2, 6);
        else 
            cv::putText(image, formatDoubleValue(theta[i], 1), cv::Point(x_max - x_min + 5 * x_add + theta[i] * scale_x + 20, y_max + 2 * y_add - s0 * scale_y + 125), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(0,69,255), 2, 6);
        cv::putText(image, formatDoubleValue(s[i], 1), cv::Point(x_max - x_min + 5 * x_add - 240, y_max + 2 * y_add - (s0+s[i]) * scale_y - 80), cv::FONT_HERSHEY_SIMPLEX, 3, cv::Scalar(0,69,255), 2, 6);

        // Add image to video:
        images.push_back(image);
        video.write(image);

    }

    for (int i = 0 ; i < 300; i++) {
        video.write(images[images.size()-1]);
    }

    // Release VideoWriter:
    video.release();
    

}


