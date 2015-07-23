#include <ros/ros.h>
#include <tf/transform_broadcaster.h>
#include <sensor_msgs/Joy.h>
#include <boost/thread.hpp>

using namespace std;

vector<string> tfs;
vector<double> xtfs;
vector<double> ytfs;

int selectedTf=0;


#define DEAD_ZONE 0.05
#define DEAD_ZONE_YAW 0.05


double threshold(double x, double x_s) {
    if(x > x_s) {
        return (x - x_s)/(1.0 - x_s);
    }
    else if(x < -x_s) {
        return -(x + x_s)/(-1.0 + x_s);
    }
    else
        return 0.0;
}



void joyCallback(const sensor_msgs::Joy::ConstPtr& msg) {
    int switchTf=msg->buttons[0];
    int selectAll=msg->buttons[5];

    if (selectAll==1) {
		for (int i=0; i<tfs.size();i++) {
    			//double newX=xtfs[i]+threshold(msg->axes[1], DEAD_ZONE);
				//double newY=ytfs[i]+threshold(msg->axes[0], DEAD_ZONE);
			double newX=xtfs[i]+threshold(msg->axes[0], DEAD_ZONE)*0.02*1.5;
			double newY=ytfs[i]+threshold(msg->axes[1], DEAD_ZONE)*0.02*1.5;

				xtfs[i]=newX;
    			ytfs[i]=newY;
    		}
    }
    else {
    	if (switchTf==1) {
    		selectedTf++;
    		if (selectedTf==tfs.size()) {
    			selectedTf=0;
    		}
           	cout<<"Switching "<<tfs[selectedTf]<<"\n";

    	}
		double newX=xtfs[selectedTf]+threshold(msg->axes[0], DEAD_ZONE)*0.02*1.5;
				double newY=ytfs[selectedTf]+threshold(msg->axes[1], DEAD_ZONE)*0.02*1.5;
    		xtfs[selectedTf]=newX;
    		ytfs[selectedTf]=newY;
    }


}

void sendPosition() {
    ros::Rate r(10);
    while (ros::ok()) {
        for (int i=0; i<tfs.size(); i++) {
         static tf::TransformBroadcaster br;
        
         tf::Transform transform;
         transform.setOrigin( tf::Vector3(xtfs[i], ytfs[i], 0.0) );
         tf::Quaternion q;
         q.setRPY(0, 0, 0);
         transform.setRotation(q);
         br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "map", tfs[i]));
        }
        r.sleep();
    }
}

int main (int argc, char** argv) {
    
    for (int i=1; i<argc; i++) {
        tfs.push_back(argv[i]);
        cout<<"input "<<argv[i]<<"\n";
        xtfs.push_back(0);
        ytfs.push_back(0);
    }

    ros::init(argc,argv,"tf_joypad");
    ros::NodeHandle n;

    boost::thread t(sendPosition);
    ros::Subscriber joy_sub= n.subscribe<sensor_msgs::Joy>("velocity_control", 10, joyCallback);

    ros::spin();
    ros::shutdown();
}



