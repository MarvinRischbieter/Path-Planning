#include <fstream>
#include <cstdlib>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

#define MAX_VEL 22.15 // 50 mph is 22.352 mps, so use lower value
#define MAX_S 6945.554 // The max s value before wrapping around the track back to 0
#define S_DIST 30  // Step size in m
#define MAX_ACC 0.15  // Maximal acceleration
#define N_POINTS 50  // Number of points of the smoothed path to generate

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for(int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x,y,map_x,map_y);
        if(dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

    int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2((map_y-y),(map_x-x));

    double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

    int prev_wp;
    prev_wp = next_wp-1;
    if(next_wp == 0)
    {
        prev_wp  = maps_x.size()-1;
    }

    double n_x = maps_x[next_wp]-maps_x[prev_wp];
    double n_y = maps_y[next_wp]-maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
    double proj_x = proj_norm*n_x;
    double proj_y = proj_norm*n_y;

    double frenet_d = distance(x_x,x_y,proj_x,proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000-maps_x[prev_wp];
    double center_y = 2000-maps_y[prev_wp];
    double centerToPos = distance(center_x,center_y,x_x,x_y);
    double centerToRef = distance(center_x,center_y,proj_x,proj_y);

    if(centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for(int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
    }

    frenet_s += distance(0,0,proj_x,proj_y);

    return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
    int prev_wp = -1;

    while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp+1)%maps_x.size();

    double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s-maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
    double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

    double perp_heading = heading-pi()/2;

    double x = seg_x + d*cos(perp_heading);
    double y = seg_y + d*sin(perp_heading);

    return {x,y};

}

// Estimate the speed of cars in a lane
double speed_of_lane(vector<vector<double>> sensor, double car_s,
                     double dt, int line, double dist){
    // Return the slowest car speed in dist interval from car_s.
    // Return MAX_VEL is no cars found
    double lane_speed = MAX_VEL;
    for (int i=0; i < sensor.size(); i++){
        double d = sensor[i][6];
        if (d < 4*(line+1) && d > 4*line){  // The car in the considered line
            //Speed of the considered car
            const double vx = sensor[i][3];
            const double vy = sensor[i][4];
            const double speed = sqrt(vx*vx+vy*vy);
            // Actual car position at the considered moment
            double s = sensor[i][5] + dt * speed;
            if (s > MAX_S){ // We have a round track!
                s -= MAX_S;
            }

            if (s > car_s && s - car_s < dist){
                // Tha car is close ahead
                if (speed < lane_speed){
                    lane_speed = speed;
                }
            }
        }
    }
    return lane_speed;
}

// Check if it is safe to change lane
bool safe_lane_change(vector<vector<double>> sensor, double car_s,
                      double dt, int line){
    bool safe_change = true;
    for (int i=0; i < sensor.size(); i++){
        double d = sensor[i][6];
        if (d < 4*(line+1) && d > 4*line){  // The car in the considered line
            const double vx = sensor[i][3];
            const double vy = sensor[i][4];
            const double speed = sqrt(vx*vx+vy*vy);
            double s = sensor[i][5] + dt * speed;
            if (s > MAX_S){
                s -= MAX_S;
            }
            double dist = s - car_s;
            if (dist<S_DIST && dist>0){
                safe_change = false;
                break;
            }
        }
    }
    return safe_change;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  double max_s = MAX_S;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }
  
  int lane = 1;  // the current line
  double cur_speed = 0; // current speed
  bool change = false;  // we are changing the lane right now

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s,
               &map_waypoints_dx, &map_waypoints_dy, &lane, &cur_speed,
               &change](uWS::WebSocket<uWS::SERVER> ws, char *data,
               size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(data);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object

            // Main car's localization Data
            double car_x = j[1]["x"];
            double car_y = j[1]["y"];
            double car_s = j[1]["s"];
            double car_d = j[1]["d"];
            double car_yaw = j[1]["yaw"];
            double car_speed = j[1]["speed"];

            // Previous path data given to the Planner
            auto previous_path_x = j[1]["previous_path_x"];
            auto previous_path_y = j[1]["previous_path_y"];
            // Previous path's end s and d values 
            double end_path_s = j[1]["end_path_s"];
            double end_path_d = j[1]["end_path_d"];

            // Sensor Fusion Data, a list of all other cars on the same side of the road.
            auto sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();
            if (prev_size > 0) {
              car_s = end_path_s;
            }

            vector<double> ptsx;
            vector<double> ptsy;

            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            
            // If previous size is almost empty, use the car as starting point
            if (prev_size < 2){
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);
                ptsx.push_back(prev_car_x);
                ptsx.push_back(car_x);
                ptsy.push_back(prev_car_y);
                ptsy.push_back(car_y);
            }
            else{
                ref_x = previous_path_x[prev_size-1];
                ref_y = previous_path_y[prev_size-1];
                double ref_x_prev = previous_path_x[prev_size-2];
                double ref_y_prev = previous_path_y[prev_size-2];
                ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
                ptsx.push_back(ref_x_prev);
                ptsx.push_back(ref_x);
                ptsy.push_back(ref_y_prev);
                ptsy.push_back(ref_y);
            }

            // Measure max speeds on all three lanes
            double lane_speeds[3];
            for (int i = 0; i < 3; i++){
                double speed_i = speed_of_lane(sensor_fusion, car_s,
                                   S_DIST/cur_speed, i, 2*S_DIST);
                lane_speeds[i] = speed_i;
            }

            // Find the fastest lane to use
            int best_lane = lane;
            for (int i = 0; i < 3; i++){
                // Consider line change only to the nearest lane
                // if it improves our speed meaningfully (at least by 0.5)
                if (lane_speeds[i] - lane_speeds[best_lane] > 0.5
                    && abs(i-lane) < 2){
                    best_lane = i;
                }
            }
            
            // Check if it is safe to change lane to the fastest one
            int next_lane = lane;
            // In case the maneuver is safe, do it!
            if (safe_lane_change(sensor_fusion,
                car_s, S_DIST/cur_speed, best_lane)){
                if (!change && best_lane != lane){
                    change = true;
                }
                next_lane = best_lane;
            }

            // Set the car speed
            double desired_speed = lane_speeds[next_lane];
            if (!change){
                // If we are not changing a lane, allow to e closer to leading car
                desired_speed = speed_of_lane(sensor_fusion, car_s,
                                   S_DIST/cur_speed, lane, S_DIST);
            }

            // Create coarse path (3 points long)
            for (int i=0; i < 3; i++){
                // Create 3 waypoints
                vector<double> wp = getXY(car_s+S_DIST*(i+1), (2 + 2*(next_lane + lane)), map_waypoints_s, map_waypoints_x, map_waypoints_y);
                ptsx.push_back(wp[0]);
                ptsy.push_back(wp[1]);
            }
            
            // Complite lane change if it was in progress
            lane = next_lane;
            if (abs(car_d - (2 + 4*next_lane)) < 0.25){
                // If the lane change was complited
                change = false;
            }

            for (int i = 0; i < ptsx.size(); i++){
                double shift_x = ptsx[i]-ref_x;
                double shift_y = ptsy[i]-ref_y;
                ptsx[i] = shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw);
                ptsy[i] = shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw);
            }

            tk::spline s;  // Define the spline
            s.set_points(ptsx, ptsy); 

            json msgJson;

            vector<double> next_x_vals;
            vector<double> next_y_vals;
            for (int i = 0; i < prev_size; i++){
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }

            double target_y = s(S_DIST);
            double target_dist = sqrt(target_y*target_y+S_DIST*S_DIST);
            double x_add_on = 0;
            // Smooth speed transition to the desired speed
            if (cur_speed < desired_speed) {
                cur_speed += MAX_ACC;
            }
            else {
                cur_speed -= MAX_ACC;
            }
            
            // Create the smoothed path with the spline
            double N = target_dist / (0.02 * cur_speed);
            double cos_y = cos(ref_yaw);
            double sin_y = sin(ref_yaw);
            for (int i = 0; i < N_POINTS-prev_size; i++){
                double x_point = x_add_on+S_DIST/N;
                double y_point = s(x_point);
                x_add_on = x_point;
                double x_ref = x_point;
                double y_ref = y_point;
                x_point = x_ref * cos_y - y_ref * sin_y + ref_x;
                y_point = x_ref * sin_y + y_ref * cos_y + ref_y;
                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }

            // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            msgJson["next_x"] = next_x_vals;
            msgJson["next_y"] = next_y_vals;

            auto msg = "42[\"control\","+ msgJson.dump()+"]";

            //this_thread::sleep_for(chrono::milliseconds(1000));
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    cout << "Connected!!!" << endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    cout << "Disconnected" << endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    cout << "Listening to port " << port << endl;
  } else {
    cerr << "Failed to listen to port" << endl;
    return -1;
  }
  h.run();
}
