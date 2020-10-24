#include<iostream>
#include<set>
#include<algorithm>
#include<vector>
#include<queue>
#include<cstdlib>
#include<cfloat>
#include<ctime>
using namespace std;

class vertex {
public:
	class edge* edge_cw;
	class edge* edge_ccw;
	double v_x;
	double v_y;
	vertex(double in_x, double in_y) { v_x = in_x; v_y = in_y; }
};

class point {//not a vertex, its position is changed by the rectangle's position
public:
	double x_c1;
	double y_c1;
	double x_c2;
	double y_c2;
	double c0;
	double c1;
	void set_value(double x_co1, double y_co1, double co0, double x_co2, double y_co2, double co1) {
		x_c1 = x_co1;
		y_c1 = y_co1;
		c0 = co0;
		x_c2 = x_co2;
		y_c2 = y_co2;
		c1 = co1;
	}
	point() {};
	point(double co0, double co1, double x_co1=0, double y_co1=0, double x_co2=0, double y_co2=0) {
		x_c1 = x_co1;
		y_c1 = y_co1;
		c0 = co0;
		x_c2 = x_co2;
		y_c2 = y_co2;
		c1 = co1;
	}
	
};

class edge {
public:
	vertex * from;
	vertex* to;
	double temp_x1;
	double temp_x2;
	double temp_y1;
	double temp_y2;
	pair<bool,double> gradient; // if gradient.first equals to true, edge is perpendicular to x-axis
	void set_coefficient(point* temp, int type, double width, double height) {//type 1: horizontal / type 2: vertical
		double x_co1;
		double y_co1;
		double c0;
		double x_co2;
		double y_co2;
		double c1;
		if (type == 1)
		{
			x_co1 = 0;
			y_co1 = (from->v_x - to->v_x) / (from->v_y - to->v_y);
			c0 = from->v_x + (from->v_y - height)*(to->v_x - from->v_x) / (from->v_y - to->v_y);
			x_co2 = 0;
			y_co2 = 1;
			c1 = height;

		}
		else {
			x_co1 = 1;
			y_co1 = 0;
			c0 = width;
			x_co2 = (from->v_y - to->v_y) / (from->v_x - to->v_x);
			y_co2 = 0;
			c1 = from->v_y + (width - from->v_x)*(from->v_y - to->v_y) / (from->v_x - to->v_x);
		}
		temp->set_value(x_co1, y_co1, c0, x_co2, y_co2, c1);
		

	}
	double x_intersect(double line_y) {// return x value which intersect with given horizontal line
		if (abs(from->v_x - to->v_x )<=DBL_EPSILON)
			return from->v_x;
		else {
			return from->v_x + (from->v_y - line_y)*(to->v_x - from->v_x) / (from->v_y - to->v_y);
		}
	}
	double y_intersect(double line_x) {
		if (abs(from->v_y - to->v_y )<= DBL_EPSILON)
			return from->v_y;
		else
			return from->v_y + (line_x - from->v_x)*(from->v_y - to->v_y) / (from->v_x - to->v_x);
	}
	pair<bool,pair<double, double>> y_check(double line_y) {//if the given horizontal line intersects with the edge, return the intersection point. if not, return (false,(0,0))
		double y_high = max(from->v_y,to->v_y);
		double y_low = min(from->v_y, to->v_y);
		if (y_high > line_y && line_y >= y_low &&(y_high!=y_low)) {
			return pair<bool, pair<double, double>>(true, pair<double, double>(x_intersect(line_y), line_y));
		}
		else {
			return pair<bool,pair<double,double>>(false, pair<double, double>(0, 0));
		}
	}
	void set_gradient() {
		if (abs(from->v_x - to->v_x)<=DBL_EPSILON)
		{
			gradient.first = true;
			gradient.second = 0;
		}
		else
		{
			gradient.first = false;
			gradient.second = (from->v_y - to->v_y) / (from->v_x - to->v_x);
		}
	}
	pair<bool, pair<double, double>> x_check(double line_x) {
		double x_high = max(from->v_x, to->v_x);
		double x_low = min(from->v_x, to->v_x);
		if (x_high > line_x && line_x >= x_low && (x_high != x_low)) {
			return pair<bool, pair<double, double>>(true, pair<double, double>(line_x,y_intersect(line_x)));
		}
		else
			return pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
	}
	int check_position(double v_x, double v_y) {
		if (abs(temp_x1 - temp_x2)<=DBL_EPSILON) {
			if (v_x < temp_x1)
				return -1;
			else if (abs(v_x - temp_x1)<=DBL_EPSILON)
				return 0;
			else
				return 1;
		}
		else {
			double value_y = (temp_y2 - temp_y1)*(v_x - temp_x1) / (temp_x2 - temp_x1) + temp_y1;
			if (value_y < v_y)
				return -1;
			else if (abs(value_y - v_y)<=DBL_EPSILON)
				return 0;
			else
				return 1;
		}
	}
};

class partial_area {
public:
	double upper_gradient;
	double upper_y_c;
	double lower_gradient;
	double lower_y_c;
	double start;
	double end;
	double area;
	partial_area(pair<double,double> u_v1,pair<double,double> u_v2, pair<double,double> l_v1,pair<double,double> l_v2, double i_start, double i_end) {
		upper_gradient = (u_v2.second-u_v1.second)/(u_v2.first-u_v1.first);
		upper_y_c = u_v1.second+(u_v2.second-u_v1.second)*u_v1.first/(u_v1.first-u_v2.first);
		lower_gradient = (l_v2.second - l_v1.second) / (l_v2.first - l_v1.first);
		lower_y_c = l_v1.second + (l_v2.second - l_v1.second)*l_v1.first / (l_v1.first - l_v2.first);
		start = i_start;
		end = i_end;
		area = (upper_gradient - lower_gradient)*(i_end*i_end - i_start * i_start)/2.0 + (upper_y_c - lower_y_c)*(i_end - i_start);
	}
	double left_partial_area(double i_start) {
		return (upper_gradient - lower_gradient)*(end*end - i_start * i_start) + (upper_y_c - lower_y_c)*(end - i_start);
	}
	double right_partial_area(double i_end) {
		return (upper_gradient - lower_gradient)*(i_end*i_end - start * start) + (upper_y_c - lower_y_c)*(i_end - start);
	}
	double find_max_area(partial_area r_area, double domain) {//compute maximum area when rectangle intersect with this and another partial_area "r_area". domain is remained width of rectangle
		double f_g = r_area.upper_gradient - r_area.lower_gradient + lower_gradient - upper_gradient;
		double f_y_c = r_area.upper_y_c - r_area.lower_y_c + lower_y_c - upper_y_c;
		double max_point = -f_y_c / f_g;
		double left_max, right_max;
		if (f_g<0 && max_point > 0 && max_point < domain && end-domain+max_point>=start && r_area.start+max_point<=r_area.end) {
				return (upper_gradient - lower_gradient)*(end*end - (end - domain + max_point)*(end - domain + max_point)) / 2 + (upper_y_c - lower_y_c)*(domain - max_point)+(r_area.upper_gradient-r_area.lower_gradient)*(2*r_area.start*max_point+max_point*max_point)/2+(r_area.upper_y_c-r_area.lower_y_c)*max_point;
		}
		else {
			double remain = (end - start)-domain;
			if (remain >= 0) {
				left_max = (upper_gradient - lower_gradient)*(end*end - (end - domain)*(end - domain)) / 2 + (upper_y_c - lower_y_c)*(domain);
			}
			else {
				left_max=area+ (r_area.upper_gradient - r_area.lower_gradient)*(2 * r_area.start*(-remain) + (-remain) * (-remain)) / 2 + (r_area.upper_y_c - r_area.lower_y_c)*(-remain);
			}
			remain = r_area.end - r_area.start - domain;
			if (remain >= 0) {
				right_max= (r_area.upper_gradient - r_area.lower_gradient)*(2 * r_area.start*domain + domain * domain) / 2 + (r_area.upper_y_c - r_area.lower_y_c)*domain;
			}
			else {
				right_max = r_area.area + (upper_gradient - lower_gradient)*(end*end - (end + remain)*(end +remain)) / 2 + (upper_y_c - lower_y_c)*(-remain);
			}
			return max(left_max, right_max);
		}
	}
};

int binary_search(vector<double> y_list, vertex* top, vertex* bot, int start, int end, double width, double height);
double find_max_area(double y1, vertex* top, vertex* bot, double width, double height);
void trans_polygon(vertex* left, vertex* right, vertex** t_top, vertex** t_bot);
void find_max_point(double x_low, double x_high, double y_low, double y_high, double width, double height, vertex* top, vertex* bot, vertex* left, vertex* right);
int edge_check(edge* e_temp, double left, double right, double top, double bot);
pair<double, double> find_intersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
bool recursive_find(double x_low, double x_high, double y_low, double y_high, double width, double height, vector<int> e_check_list, edge* e_list[], int d_num, vector<pair<double, double>> domain);
void cal_triangle(double *c_x_2, double *c_x, double* c_xy, double* c_y, double* c_y_2, double* con, point* a, point* b, point* c,double test_x, double test_y);
void product_point(double* c_x_2, double* c_x, double* c_xy, double* c_y, double* c_y_2, double* con, point* a, point* b);
bool check_point(double p_x, double p_y, vector<pair<double, double>>domain);
void find_max_area_side(double x_2, double xy, double y_2, double x_1, double y_1, double constant, vector<pair<double, double>> domain);
pair<double, double> max_pos;
double global_max = 0;
bool global_check = false;
int main() {
	

	int num_v;
	
	vertex* top_vertex;
	vertex* bot_vertex;
	vertex* left_vertex;
	vertex* right_vertex;
	double left_side,right_side,top_side,bottom_side;
	vertex* v_temp;
	vertex* start_node;
	vertex* t_bot=NULL;
	vertex* t_top=NULL;
	double width, height;
	cout << "width and heigth:";
	cin >> width >> height;
	double temp_x, temp_y;
	cin >> num_v;
	cin >> temp_x >> temp_y;//order must be counter clock wise
	left_side = right_side = temp_x;
	top_side = bottom_side = temp_y;

	v_temp = new vertex(temp_x, temp_y);
	start_node = v_temp;
	top_vertex = v_temp;
	bot_vertex = v_temp;
	left_vertex = v_temp;
	right_vertex = v_temp;

	//get user input
	for (int i = 1; i < num_v; i++) {
		cin >> temp_x >> temp_y;

		left_side = min(left_side, temp_x);
		right_side = max(right_side, temp_x);
		top_side = max(top_side, temp_y);
		bottom_side = min(bottom_side, temp_y);

		v_temp->edge_ccw = new edge();
		v_temp->edge_ccw->from = v_temp;
		v_temp->edge_ccw->to = new vertex(temp_x, temp_y);
		v_temp->edge_ccw->set_gradient();
		if (top_side == temp_y)
			top_vertex = v_temp->edge_ccw->to;
		if (bottom_side == temp_y)
			bot_vertex = v_temp->edge_ccw->to;
		if (left_side == temp_x)
			left_vertex = v_temp->edge_ccw->to;
		if (right_side == temp_x)
			right_vertex = v_temp->edge_ccw->to;
		v_temp->edge_ccw->to->edge_cw = v_temp->edge_ccw;
		v_temp = v_temp->edge_ccw->to;
		if (i ==num_v-1) {
			v_temp->edge_ccw = new edge();
			v_temp->edge_ccw->from = v_temp;
			v_temp->edge_ccw->to = start_node;
			v_temp->edge_ccw->set_gradient();
			start_node->edge_cw = v_temp->edge_ccw;
		}
	
	}
	clock_t start = clock();
	// create x_list and y_list
	vertex* l_temp = bot_vertex->edge_cw->from;
	vertex* r_temp = bot_vertex->edge_ccw->to;
	vertex* u_temp = left_vertex->edge_cw->from;
	vertex* d_temp = left_vertex->edge_ccw->to;
	vector<double> y_list;
	vector<double> x_list;
	y_list.push_back(bot_vertex->v_y);
	x_list.push_back(left_vertex->v_x);
	while (u_temp != right_vertex || d_temp != right_vertex) {
		if (u_temp == right_vertex) {
			while (d_temp != right_vertex) {
				x_list.push_back(d_temp->v_x);
				d_temp = d_temp->edge_ccw->to;
			}
			break;
		}
		else if (d_temp == right_vertex) {
			while (u_temp != right_vertex) {
				x_list.push_back(u_temp->v_x);
				u_temp = u_temp->edge_cw->from;
			}
			break;
		}
		if (u_temp->v_x > d_temp->v_x) {
			x_list.push_back(d_temp->v_x);
			d_temp = d_temp->edge_ccw->to;
		}
		else if (u_temp->v_x < d_temp->v_x) {
			x_list.push_back(u_temp->v_x);
			u_temp = u_temp->edge_cw->from;
		}
		else {
			x_list.push_back(d_temp->v_x);
			u_temp = u_temp->edge_cw->from;
			d_temp = d_temp->edge_ccw->to;
		}
	}

	while (l_temp != top_vertex || r_temp != top_vertex) {
		if (l_temp == top_vertex) {
			while (r_temp != top_vertex) {
				y_list.push_back(r_temp->v_y);
				r_temp = r_temp->edge_ccw->to;
			}
			break;
		}
		else if (r_temp == top_vertex) {
			while (l_temp != top_vertex) {
				y_list.push_back(l_temp->v_y);
				l_temp = l_temp->edge_cw->from;
			}
			break;
		}
		if (l_temp->v_y > r_temp->v_y) {
			y_list.push_back(r_temp->v_y);
			r_temp = r_temp->edge_ccw->to;
		}
		else if (l_temp->v_y < r_temp->v_y) {
			y_list.push_back(l_temp->v_y);
			l_temp = l_temp->edge_cw->from;
		}
		else {
			y_list.push_back(r_temp->v_y);
			l_temp = l_temp->edge_cw->from;
			r_temp = r_temp->edge_ccw->to;
		}
	}
	y_list.push_back(top_vertex->v_y);
	x_list.push_back(right_vertex->v_x);
	vector<double> y_list_prime;
	vector<double> x_list_prime;
	vector<double> y_list_res;
	vector<double> x_list_res;
	vector<double>::iterator it,it2;
	for (it = x_list.begin(); it != x_list.end(); it++)
		x_list_prime.push_back(*it - width);
	for (it = y_list.begin(); it != y_list.end(); it++) {
		y_list_prime.push_back(*it - height);
	}
	it = x_list.begin();
	it2 = x_list_prime.begin();
	while (it != x_list.end() || it2 != x_list_prime.end()) {
		if (it == x_list.end()) {
			while (it2 != x_list_prime.end()) {
				x_list_res.push_back(*it2);
				it2++;
			}
			break;
		}
		else if (it2 == x_list_prime.end()) {
			while (it != x_list.end()) {
				x_list_res.push_back(*it);
				it++;
			}
			break;
		}

		if (*it < *it2)
		{
			x_list_res.push_back(*it);
			it++;
		}
		else if (*it > *it2)
		{
			x_list_res.push_back(*it2);
			it2++;
		}
		else {
			x_list_res.push_back(*it);
			it++;
			it2++;
		}
	}

	it = y_list.begin();
	it2 = y_list_prime.begin();
	while (it != y_list.end() || it2 != y_list_prime.end()) {
		if (it == y_list.end()) {
			while (it2 != y_list_prime.end()) {
				y_list_res.push_back(*it2);
				it2++;
			}
			break;
		}
		else if (it2 == y_list_prime.end()) {
			while (it != y_list.end()) {
				y_list_res.push_back(*it);
				it++;
			}
			break;
		}

		if (*it < *it2)
		{
			y_list_res.push_back(*it);
			it++;
		}
		else if (*it > *it2)
		{
			y_list_res.push_back(*it2);
			it2++;
		}
		else {
			y_list_res.push_back(*it);
			it++;
			it2++;
		}
	}

	int y_idx;
	int x_idx;

	trans_polygon(left_vertex, right_vertex, &t_top, &t_bot);

	//find y_list_res's index which has maximum area on horizontal line with y value
	y_idx = binary_search(y_list_res, top_vertex, bot_vertex, 0, y_list_res.size()-1, width, height);
	
	x_idx = binary_search(x_list_res, t_top, t_bot, 0, x_list_res.size()-1, height, width);
	

	if (x_idx != 0) {
		if(y_idx!=0)
			find_max_point(x_list_res[x_idx - 1], x_list_res[x_idx], y_list_res[y_idx - 1], y_list_res[y_idx], width, height, top_vertex, bot_vertex, left_vertex, right_vertex);
		if(y_idx!=y_list_res.size()-1)
			find_max_point(x_list_res[x_idx - 1], x_list_res[x_idx], y_list_res[y_idx], y_list_res[y_idx + 1], width, height, top_vertex, bot_vertex, left_vertex, right_vertex);
	}
	
	if(x_idx != x_list_res.size()-1 && y_idx!=0)
		find_max_point(x_list_res[x_idx], x_list_res[x_idx + 1], y_list_res[y_idx - 1], y_list_res[y_idx], width, height, top_vertex, bot_vertex, left_vertex, right_vertex);
	if(x_idx!=x_list_res.size()-1 && y_idx!=y_list_res.size()-1)
		find_max_point(x_list_res[x_idx], x_list_res[x_idx + 1], y_list_res[y_idx], y_list_res[y_idx + 1], width, height, top_vertex, bot_vertex, left_vertex, right_vertex);
	
	
	cout.precision(8);
	cout << "answer: (x,y) = (" << max_pos.first << " , " << max_pos.second << " ) /// area: " << global_max << endl;
	
	clock_t end = clock();
	
	cout <<  "\nTime: " << (double)(end - start)/CLOCKS_PER_SEC  << endl;
	return 0;
}
void find_max_point(double x_low,double x_high, double y_low, double y_high, double width, double height, vertex* top, vertex* bot,vertex* left, vertex* right) {
	edge* e_list[8] = { NULL };//the list including the edges that possibly intersect with rectangle. from the left bottom vertex, ccw
	vector<int> e_check_list;// 0: not yet decided, 1: include, 2: not include 3: not intersect
	for (int i = 0; i < 8; i++) {
		e_check_list.push_back(0);
	}
	vertex* temp = bot;
	vector<pair<double, double>> domain;
	edge* e_temp = bot->edge_ccw;
	int phase = 1;
	int d_num=0;
	pair<bool, pair<double, double>> p_temp;
	
	 do{
		if (!edge_check(e_temp, x_low, x_high + width, y_high + height, y_low))
		{
			temp = temp->edge_ccw->to;
			e_temp = temp->edge_ccw;
			if (temp == right)
				phase = 2;
			else if (temp == top)
				phase = 3;
			else if (temp == left)
				phase = 4;
			continue;
		}
		else {
			p_temp = e_temp->y_check((y_low + y_high) / 2);
			if (p_temp.second.first >= x_high + width || p_temp.second.first <= x_low)
				p_temp.first = false;
			if (p_temp.first) {
				if (e_temp->gradient.first) {
					if (e_temp->from->v_x == left->v_x)
					{
						e_list[0] = e_temp;
						e_list[5] = e_temp;
					}
					else {
						e_list[1] = e_temp;
						e_list[4] = e_temp;
					}
				}
				else {
					if (phase >=3)
						e_list[0] = e_temp;
					else
						e_list[1] = e_temp;
				}

			}
			p_temp = e_temp->y_check((y_low + y_high) / 2 + height);
			if (p_temp.second.first >= x_high + width || p_temp.second.first <= x_low)
				p_temp.first = false;
			if (p_temp.first){
				if (e_temp->gradient.first != true) {
					if (phase >= 3)
						e_list[5] = e_temp;
					else
						e_list[4] = e_temp;
				}
			}

			p_temp = e_temp->x_check((x_low + x_high) / 2);
			if (p_temp.second.second >= y_high + height || p_temp.second.second <= y_low)
				p_temp.first = false;
			if (p_temp.first) {
				if (abs(e_temp->gradient.second) <= FLT_EPSILON) {
					if (abs(e_temp->from->v_y - top->v_y)<=FLT_EPSILON) {
						e_list[3] = e_temp;
						e_list[6] = e_temp;
					}
					else {
						e_list[2] = e_temp;
						e_list[7] = e_temp;
					}
				}
				else {
					if (phase == 2 || phase == 3)
						e_list[6] = e_temp;
					else
						e_list[7] = e_temp;
				}
			}
			p_temp = e_temp->x_check((x_low + x_high) / 2 + width);
			if (p_temp.second.second >= y_high + height || p_temp.second.second <= y_low)
				p_temp.first = false;
			if (p_temp.first) {
				if (abs(e_temp->gradient.second) >= DBL_EPSILON) {
					if (phase == 2 || phase == 3)
						e_list[3] = e_temp;
					else
						e_list[2] = e_temp;
				}
			}
		}
		temp = temp->edge_ccw->to;
		e_temp = temp->edge_ccw;
		if (temp == right)
			phase = 2;
		else if (temp == top)
			phase = 3;
		else if (temp == left)
			phase = 4;
		
	 } while (temp != bot);
	 for (int i = 0; i < 8; i++) {
		 if (e_list[i] == NULL)
		 {
			 e_check_list[i] = 3;
			 d_num++;
		 }
	 }
	 domain.push_back(pair<double, double>(x_low, y_low));
	 domain.push_back(pair<double, double>(x_high, y_low));
	 domain.push_back(pair<double, double>(x_high, y_high));
	 domain.push_back(pair<double, double>(x_low, y_high));
	 recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
}

void cal_triangle(double *c_x_2,double *c_x, double* c_xy, double* c_y,double* c_y_2, double* con, point* a, point* b, point* c, double test_x, double test_y) {
	double temp1 = 0;
	double temp2 = 0;
	double temp3 = 0;
	double temp4 = 0;
	double temp5 = 0;
	double temp6 = 0;

	product_point(&temp1, &temp2, &temp3, &temp4, &temp5, &temp6, a, b);
	product_point(&temp1, &temp2, &temp3, &temp4, &temp5, &temp6, b, c);
	product_point(&temp1, &temp2, &temp3, &temp4, &temp5, &temp6, c, a);

	temp1 /= 2.0;
	temp2 /= 2.0;
	temp3 /= 2.0;
	temp4 /= 2.0;
	temp5 /= 2.0;
	temp6 /= 2.0;

	/*if (temp1*test_x*test_x + temp3 * test_x* test_y + temp5 * test_y*test_y + temp2 * test_x + temp4 * test_y + temp6 < -DBL_EPSILON )
	{
		temp1 = -temp1;
		temp2 = -(temp2);
		temp3 = -(temp3);
		temp4 = -(temp4);
		temp5 = -(temp5);
		temp6 = -(temp6);
	}*/

	(*c_x_2) += temp1;
	(*c_x) += temp2;
	(*c_xy) += temp3;
	(*c_y) += temp4;
	(*c_y_2) += temp5;
	(*con) += temp6;
	
}

void product_point(double* c_x_2, double* c_x, double* c_xy, double* c_y, double* c_y_2,double* con,point* a, point* b) {
	(*c_x_2) += a->x_c1*b->x_c2 - a->x_c2*b->x_c1;
	(*c_xy) += (a->x_c1*b->y_c2) + (a->y_c1*b->x_c2) - a->x_c2*b->y_c1-a->y_c2*b->x_c1;
	(*c_x) += (a->x_c1*b->c1) + (a->c0*b->x_c2) - a->x_c2*b->c0 - a->c1*b->x_c1;
	(*c_y) += (a->y_c1*b->c1) + (a->c0*b->y_c2)-a->y_c2*b->c0 - a->c1* b->y_c1;
	(*con) += a->c0*b->c1 - a->c1*b->c0;
	(*c_y_2) += a->y_c1*b->y_c2 - a->y_c2*b->y_c1;
}

bool check_point(double p_x, double p_y, vector<pair<double, double>>domain) {
	pair<double, double> temp;
	pair<double, double> temp2;
	int check = 0; 
	int i = 0;
	for (i = 0; i < domain.size()-1; i++) {
		temp = pair<double, double>(domain[i+1].first-domain[i].first, domain[i+1].second-domain[i].second);
		temp2 = pair<double, double>(p_x - domain[i].first, p_y - domain[i].second);
		if (temp.first*temp2.second - temp2.first*temp.second < 0)
			return false;

	}
	temp = pair<double, double>(domain[0].first - domain[i].first, domain[0].second - domain[i].second);
	temp2 = pair<double, double>(p_x - domain[i].first, p_y - domain[i].second);
	if (temp.first*temp2.second - temp2.first*temp.second < 0)
		return false;
	return true;
}

bool recursive_find(double x_low, double x_high, double y_low, double y_high, double width, double height, vector<int> e_check_list, edge* e_list[], int d_num, vector<pair<double, double>> domain) {
	int i = 0;
	int k = 0;
	if (d_num == 8) {
	
		i = 0;
		k = -1;
		int order = 0;
		int check = 0;
		vector<int> inter_vector;
		while (i < 8) {
			if (e_check_list[i] == 1)
			{
				if (order == 1 && i % 2 == 0 || order == -1 && i % 2 == 1)
					return false;
				if (i % 2 == 0)
					order = 1;
				else
					order = -1;
				inter_vector.push_back(i);
				check++;
				
			}
			i++;
		}
		if (check % 2 == 1)
			return false;
		// this is the case that rectangle can be included in convex polygon
		if (check == 0) {
			cout << "answer1: (x,y) = (" << domain[0].first << " , " << domain[0].second << " ) /// area: " << width * height << endl;
			return true;
		}

		
		vector<point*> point_v;
		point* point_temp;
		edge* edge_temp = e_list[inter_vector[0]];
		double c_x=0;
		double c_y = 0;
		double c_xy = 0;
		double c_x_2 = 0;
		double c_y_2 = 0;
		double constant = 0;
		double sol_x;
		double sol_y;
		double temp_max_value;
		int t;
		for (int j = 0; j < inter_vector.size()-1; j++) {
			point_temp = new point;
			if (inter_vector[j] < 2 )
				e_list[inter_vector[j]]->set_coefficient(point_temp,1, 0, 0);
			else if(inter_vector[j]>=6)
				e_list[inter_vector[j]]->set_coefficient(point_temp,2, 0, 0);
			else if(inter_vector[j]>=4)
				e_list[inter_vector[j]]->set_coefficient(point_temp,1, 0, height);
			else
				e_list[inter_vector[j]]->set_coefficient(point_temp,2, width, 0);
			point_v.push_back(point_temp);
			t = inter_vector[j];
			if (inter_vector[j] % 2 == 0)
			{
				t++;
				while (t != inter_vector[j + 1])
				{
					t++;
					if (t > 7)
						break;
					if (t == 2)
					{
						point_temp = new point;
						point_temp->set_value(1, 0, width, 0, 1, 0);
						point_v.push_back(point_temp);
					}
					else if (t == 4) {
						point_temp = new point;
						point_temp->set_value(1, 0, width, 0, 1, height);
						point_v.push_back(point_temp);
					}
					else if (t == 6) {
						point_temp = new point;
						point_temp->set_value(1, 0, 0, 0, 1, height);
						point_v.push_back(point_temp);
					}
					
					
				}
			}
			else {
				edge_temp = e_list[inter_vector[j]];
				while (edge_temp != e_list[inter_vector[j + 1]]) {
					point_temp = new point;
					point_temp->set_value(0, 0, edge_temp->to->v_x, 0, 0, edge_temp->to->v_y);
					point_v.push_back(point_temp);
					edge_temp = edge_temp->to->edge_ccw;
				}
			}
		}
		point_temp = new point;
		if (inter_vector[inter_vector.size()-1] < 2)
			e_list[inter_vector[inter_vector.size()-1]]->set_coefficient(point_temp, 1, 0, 0);
		else if (inter_vector[inter_vector.size() - 1] >= 6)
			e_list[inter_vector[inter_vector.size() - 1]]->set_coefficient(point_temp, 2, 0, 0);
		else if (inter_vector[inter_vector.size() - 1] >= 4)
			e_list[inter_vector[inter_vector.size() - 1]]->set_coefficient(point_temp, 1, 0, height);
		else
			e_list[inter_vector[inter_vector.size() - 1]]->set_coefficient(point_temp, 2, width, 0);
		point_v.push_back(point_temp);
		if (inter_vector[inter_vector.size() - 1] % 2 == 0)
		{
			t = inter_vector[inter_vector.size() - 1];
			t++;
			while (t < 8)
			{
				if (t == 2)
				{
					point_temp = new point;
					point_temp->set_value(1, 0, width, 0, 1, 0);
					point_v.push_back(point_temp);
				}
				else if (t == 4) {
					point_temp = new point;
					point_temp->set_value(1, 0, width, 0, 1, height);
					point_v.push_back(point_temp);
				}
				else if (t == 6) {
					point_temp = new point;
					point_temp->set_value(1, 0, 0, 0, 1, height);
					point_v.push_back(point_temp);
				}
				
				
				t++;
			}
			t = 0;
			while (t < inter_vector[0]) {
				if (t == 0) {
					point_temp = new point;
					point_temp->set_value(1, 0, 0, 0, 1, 0);
					point_v.push_back(point_temp);
				}
				else if (t == 2)
				{
					point_temp = new point;
					point_temp->set_value(1, 0, width, 0, 1, 0);
					point_v.push_back(point_temp);
				}
				else if (t == 4) {
					point_temp = new point;
					point_temp->set_value(1, 0, width, 0, 1, height);
					point_v.push_back(point_temp);
				}
				else if (t == 6) {
					point_temp = new point;
					point_temp->set_value(1, 0, 0, 0, 1, height);
					point_v.push_back(point_temp);
				}
			

				t++;
			}
		}
		else {
			edge_temp = e_list[inter_vector[inter_vector.size()-1]];
			while (edge_temp != e_list[inter_vector[0]]) {
				point_temp = new point;
				point_temp->set_value(0, 0, edge_temp->to->v_x, 0, 0, edge_temp->to->v_y);
				point_v.push_back(point_temp);
				edge_temp = edge_temp->to->edge_ccw;
			}
		}
		for (int j = 1; j < point_v.size() - 1; j++) {
			cal_triangle(&c_x_2, &c_x, &c_xy, &c_y, &c_y_2, &constant, point_v[0], point_v[j], point_v[j + 1],domain[0].first,domain[0].second);
		}
		/* test code
		if (check_point(0.876685, 1.28169, domain))
		{
			double xxx = 0.893972;
			double yyy = 1.56952;
			cout << "ysy" << endl;
			cout << c_x_2 * xxx*xxx + c_y_2 * yyy*yyy + c_xy * xxx*yyy + c_x * xxx + c_y * yyy + constant;
			cout << "check" << endl;
		}
		*/

		if (c_x_2 >= 0) {
			find_max_area_side(c_x_2, c_xy, c_y_2, c_x, c_y, constant, domain);
			return false;
		}
		else {
			if (4 * c_x_2*c_y_2 - c_xy*c_xy > 0 && c_x_2<0) {
				if (abs(c_xy) <= DBL_EPSILON) {
					sol_x = -c_x / (2 * c_x_2);
					sol_y = -c_y / (2 * c_y_2);
				}
				else {
					sol_x = (2 * c_y_2*c_x / c_xy - c_y) / (c_xy - 4 * c_x_2*c_y_2 / c_xy);
					sol_y = (2 * c_x_2*sol_x + c_x) / (-c_xy);
				}

				if (!check_point(sol_x, sol_y, domain))
				{
					find_max_area_side(c_x_2, c_xy, c_y_2, c_x, c_y, constant, domain);
					return false;
				}

				global_check = true;
				temp_max_value = c_x_2 * sol_x*sol_x + c_xy * sol_x*sol_y + c_y_2 * sol_y*sol_y + c_x * sol_x + c_y * sol_y + constant;
				
				if (global_max < temp_max_value) {
					global_max = temp_max_value;
					max_pos.first = sol_x;
					max_pos.second = sol_y;
				}
				
			}
			else {
				find_max_area_side(c_x_2, c_xy, c_y_2, c_x, c_y, constant, domain);
				return false;
			}
		}

		return true;
	}

	///////////////////////////if d_num<8
	vector<pair<double, double>> v_list1;
	vector<pair<double, double>> v_list2;
	vector<pair<double, double>> v_list3;
	vector<pair<double, double>> v_list4;
	edge e_temp;
	edge e_temp2;
	pair < double, double > i_temp;
	int pos_temp;
	int pos_prev = 0;
	int pos_first;
	double gradient;
	bool res1, res2, res3;
	pair<bool, pair<double, double>> p_temp;
	while (e_check_list[i] != 0)
	{
		i++;
		if (i == 8)
			break;
	}
	if (i == 0 || i == 1) {
		e_temp.temp_x1 = e_list[i]->from->v_x;
		e_temp.temp_y1 = e_list[i]->from->v_y;
		e_temp.temp_x2 = e_list[i]->to->v_x;
		e_temp.temp_y2 = e_list[i]->to->v_y;
		e_temp2.temp_x1 = e_list[i]->from->v_x;
		e_temp2.temp_y1 = e_list[i]->from->v_y;
		e_temp2.temp_x2 = e_list[i]->to->v_x;
		e_temp2.temp_y2 = e_list[i]->to->v_y;
		gradient = (e_temp.temp_y2 - e_temp.temp_y1) / (e_temp.temp_x2 - e_temp.temp_x1);
		if (gradient > 0.0) {
			e_temp.temp_x1 -= width;
			e_temp.temp_x2 -= width;
		}
		else {
			e_temp2.temp_x1 -= width;
			e_temp2.temp_x2 -= width;
		}
	}
	else if (i == 2 || i == 3) {
		e_temp.temp_x1 = e_list[i]->from->v_x - width;
		e_temp.temp_y1 = e_list[i]->from->v_y;
		e_temp.temp_x2 = e_list[i]->to->v_x - width;
		e_temp.temp_y2 = e_list[i]->to->v_y;
		e_temp2.temp_x1 = e_list[i]->from->v_x - width;
		e_temp2.temp_y1 = e_list[i]->from->v_y - height;
		e_temp2.temp_x2 = e_list[i]->to->v_x - width;
		e_temp2.temp_y2 = e_list[i]->to->v_y - height;

	}
	else if (i == 4 || i == 5) {
		e_temp.temp_x1 = e_list[i]->from->v_x;
		e_temp.temp_y1 = e_list[i]->from->v_y - height;
		e_temp.temp_x2 = e_list[i]->to->v_x;
		e_temp.temp_y2 = e_list[i]->to->v_y - height;
		e_temp2.temp_x1 = e_list[i]->from->v_x;
		e_temp2.temp_y1 = e_list[i]->from->v_y - height;
		e_temp2.temp_x2 = e_list[i]->to->v_x;
		e_temp2.temp_y2 = e_list[i]->to->v_y - height;
		gradient = (e_temp.temp_y2 - e_temp.temp_y1) / (e_temp.temp_x2 - e_temp.temp_x1);
		if (gradient > 0) {
			e_temp.temp_x1 -= width;
			e_temp.temp_x2 -= width;
		}
		else {
			e_temp2.temp_x1 -= width;
			e_temp2.temp_x2 -= width;
		}
	}
	else {
		e_temp.temp_x1 = e_list[i]->from->v_x;
		e_temp.temp_y1 = e_list[i]->from->v_y;
		e_temp.temp_x2 = e_list[i]->to->v_x;
		e_temp.temp_y2 = e_list[i]->to->v_y;
		e_temp2.temp_x1 = e_list[i]->from->v_x;
		e_temp2.temp_y1 = e_list[i]->from->v_y - height;
		e_temp2.temp_x2 = e_list[i]->to->v_x;
		e_temp2.temp_y2 = e_list[i]->to->v_y - height;
	}
	
	for (int k = 0; k < domain.size(); k++) {

		pos_temp = e_temp.check_position(domain[k].first, domain[k].second);// if vertex locate on the left or under, return -1 , on the edge: return 0 , upper or right: return 1
		if (k == 0)
		{
			pos_first = pos_temp;
		}

		if (pos_temp == 1 && pos_prev == -1 || pos_temp == -1 && pos_prev == 1) {
			i_temp = find_intersect(domain[k - 1].first, domain[k - 1].second, domain[k].first, domain[k].second, e_temp.temp_x1, e_temp.temp_y1, e_temp.temp_x2, e_temp.temp_y2);
			
			
			v_list1.push_back(pair<double, double>(i_temp.first, i_temp.second));
			v_list2.push_back(pair<double, double>(i_temp.first, i_temp.second));
		}
		if (pos_temp == -1)
			v_list1.push_back(pair<double, double>(domain[k].first, domain[k].second));
		else if (pos_temp == 0) {
			v_list1.push_back(pair<double, double>(domain[k].first, domain[k].second));
			v_list2.push_back(pair<double, double>(domain[k].first, domain[k].second));
		}
		else
			v_list2.push_back(pair<double, double>(domain[k].first, domain[k].second));

		pos_prev = pos_temp;
	}
	
	if (pos_first == 1 && pos_prev == -1 || pos_first == -1 && pos_prev == 1) {
		i_temp = find_intersect(domain[domain.size() - 1].first, domain[domain.size() - 1].second, domain[0].first, domain[0].second, e_temp.temp_x1, e_temp.temp_y1, e_temp.temp_x2, e_temp.temp_y2);
		v_list1.push_back(pair<double, double>(i_temp.first, i_temp.second));
		v_list2.push_back(pair<double, double>(i_temp.first, i_temp.second));
	}

	pos_prev = 0;

	for (int k = 0; k < v_list2.size(); k++) {

		pos_temp = e_temp2.check_position(v_list2[k].first, v_list2[k].second);// if vertex locate on the left or under, return -1 , on the edge: return 0 , upper or right: return 1
		if (k == 0)
		{
			pos_first = pos_temp;
		}
		if (pos_temp == 1 && pos_prev == -1 || pos_temp == -1 && pos_prev == 1) {
			i_temp = find_intersect(v_list2[k - 1].first, v_list2[k - 1].second, v_list2[k].first, v_list2[k].second, e_temp2.temp_x1, e_temp2.temp_y1, e_temp2.temp_x2, e_temp2.temp_y2);
			v_list3.push_back(pair<double, double>(i_temp.first, i_temp.second));
			v_list4.push_back(pair<double, double>(i_temp.first, i_temp.second));
		}
		if (pos_temp == -1)
			v_list3.push_back(pair<double, double>(v_list2[k].first, v_list2[k].second));
		else if (pos_temp == 0) {
			v_list3.push_back(pair<double, double>(v_list2[k].first, v_list2[k].second));
			v_list4.push_back(pair<double, double>(v_list2[k].first, v_list2[k].second));
		}
		else
			v_list4.push_back(pair<double, double>(v_list2[k].first, v_list2[k].second));

		pos_prev = pos_temp;
	}
	if (pos_first == 1 && pos_prev == -1 || pos_first == -1 && pos_prev == 1) {
		i_temp = find_intersect(v_list2[v_list2.size() - 1].first, v_list2[v_list2.size() - 1].second, v_list2[0].first, v_list2[0].second, e_temp2.temp_x1, e_temp2.temp_y1, e_temp2.temp_x2, e_temp2.temp_y2);
		v_list3.push_back(pair<double, double>(i_temp.first, i_temp.second));
		v_list4.push_back(pair<double, double>(i_temp.first, i_temp.second));
	}

	res1 = false;
	res2 = false;
	res3 = false;
	if (v_list1.size() >= 3) {
		e_check_list[i] = 2;
		res1 = recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num + 1, v_list1);
	}
	
	if (v_list3.size() >= 3) {
		e_check_list[i] = 1;
		res2 = recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num + 1, v_list3);
	}
	
	if (v_list4.size() >= 3) {
		e_check_list[i] = 2;
		res3 = recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num + 1, v_list4);
	}
	

	return false;
}

pair<double, double> find_intersect(double x1, double y1, double x2, double y2,double x3, double y3, double x4, double y4) {
	double a, b, c, d;
	if (abs(x1 - x2)<=DBL_EPSILON) {
		c = (y4 - y3) / (x4 - x3);
		return pair<double, double>(x1, c*(x1 - x3) + y3);
	}
	else if (abs(x3 - x4)<=DBL_EPSILON) {
		a = (y2 - y1) / (x2 - x1);
		return pair<double, double>(x3, a*(x3 - x1) + y1);
	}
	else {
		a = (y2 - y1) / (x2 - x1);
		b = y1 - a * x1;
		c = (y4 - y3) / (x4 - x3);
		d = y3 - c * x3;
		return pair<double, double>((d - b) / (a - c), a*(d - b) / (a - c) + b);
	}
}

int edge_check(edge* e_temp,  double left, double right, double top, double bot) {
	if ((e_temp->from->v_x <= left && e_temp->to->v_x <= left) || (e_temp->from->v_x >= right && e_temp->to->v_x >= right) || (e_temp->from->v_y <= bot && e_temp->to->v_y <= bot) || (e_temp->from->v_y >= top && e_temp->to->v_y >= top))
		return 0;
	else
		return 1;
}

double find_max_area(double y1, vertex* top, vertex* bot, double width, double height) {
	if (y1 >= top->v_y)
		return 0;

	pair<bool, pair<double, double>> bot_left, bot_right, top_left, top_right;
	bot_left = pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
	bot_right = pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
	top_left = pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
	top_right = pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
	vertex* temp = bot;
	edge* e_temp;
	
	vector<pair<double, double>> left_vertex;
	vector<pair<double, double>> right_vertex;
	vector<pair<double, double>> upper_vertex;
	vector<pair<double, double>> lower_vertex;


	////// find left vertex
	while (1) {
		e_temp = temp->edge_cw;

		
		
		if (bot_left.first == false) {
			bot_left = e_temp->y_check(y1);
			if (bot_left.first == true)
				left_vertex.push_back(pair<double, double>(bot_left.second.first, bot_left.second.second));
		}
		if (top_left.first == false) {
			top_left = e_temp->y_check(y1 + height);
			if (top_left.first == true)
			{
				if(abs(bot_left.second.second-temp->v_y)>DBL_EPSILON)
					left_vertex.push_back(pair<double, double>(temp->v_x, temp->v_y));
				if(abs(temp->v_y-top_left.second.second)>DBL_EPSILON)
					left_vertex.push_back(pair<double, double>(top_left.second.first, top_left.second.second));
			}
			else if(top_left.first == false && (bot_left.first == true || y1 <= bot->v_y+DBL_EPSILON)&&(temp->v_y>bot_left.second.second)&&(abs(temp->v_y-bot_left.second.second)>DBL_EPSILON)) {
				left_vertex.push_back(pair<double, double>(temp->v_x, temp->v_y));
			}
		}

		temp = e_temp->from;
		if (e_temp->to == top || (top_left.first&&bot_left.first))
			break;
	}
	////////////
	temp = bot;
	///////// find right vertex
	while (1) {
		e_temp = temp->edge_ccw;

		if (bot_right.first == false) {
			bot_right = e_temp->y_check(y1);
			if (bot_right.first == true)
				right_vertex.push_back(pair<double, double>(bot_right.second.first, bot_right.second.second));
		}
		if (top_right.first == false) {
			top_right = e_temp->y_check(y1 + height);
			if (top_right.first == true)
			{
				if(abs(bot_right.second.second-temp->v_y)>DBL_EPSILON)
					right_vertex.push_back(pair<double, double>(temp->v_x, temp->v_y));
				if (abs(temp->v_y - top_right.second.second)>DBL_EPSILON)
					right_vertex.push_back(pair<double, double>(top_right.second.first, top_right.second.second));
			}
			else if (top_right.first == false && (bot_right.first == true || y1 <= bot->v_y+DBL_EPSILON)&&(temp->v_y>bot_right.second.second)&&(abs(temp->v_y-bot_right.second.second)>DBL_EPSILON)) {
				right_vertex.push_back(pair<double, double>(temp->v_x, temp->v_y));
			}
		}
		temp = e_temp->to;
		if (e_temp->from == top||(top_right.first&&bot_right.first))
			break;
	}
	//////////////


	//////find left most vertex and right most vertex
	vector<pair<double, double>>::iterator it;
	double temp_min, temp_max;
	int i = 0;
	int min_idx=0, max_idx=0;
	temp_min = left_vertex.begin()->first;
	temp_max = right_vertex.begin()->first;
	for (it = left_vertex.begin(); it != left_vertex.end(); it++) {
		if (temp_min > it->first)
		{
			temp_min = it->first;
			min_idx = i;
		}
		i++;
	}
	i = 0;
	for (it = right_vertex.begin(); it != right_vertex.end(); it++) {
		if (temp_max < it->first) {
			temp_max = it->first;
			max_idx = i;
		}
		i++;
	}
	////////////

	//////seperate upper and lower vertex
	for (i = min_idx; i >= 0; i--) {
		lower_vertex.push_back(left_vertex[i]);
	}
	for (i = min_idx; i < left_vertex.size(); i++) {
		upper_vertex.push_back(left_vertex[i]);
	}
	if (abs(lower_vertex[lower_vertex.size() - 1].first - right_vertex[0].first)<=DBL_EPSILON && abs(lower_vertex[lower_vertex.size() - 1].second - right_vertex[0].second)<=DBL_EPSILON)
		i = 1;
	else
		i = 0;
	for (; i <= max_idx; i++) {
		
		lower_vertex.push_back(right_vertex[i]);
	}
	if (abs(upper_vertex[upper_vertex.size() - 1].first - right_vertex[right_vertex.size() - 1].first)<=DBL_EPSILON && abs(upper_vertex[upper_vertex.size() - 1].second - right_vertex[right_vertex.size() - 1].second)<=DBL_EPSILON)
		i = right_vertex.size() - 2;
	else
		i = right_vertex.size() - 1;
	for (; i >= max_idx; i--) {
		upper_vertex.push_back(right_vertex[i]);
	}

	////////////////
	vector<double> x_list1;
	vector<double> x_list2;
	vector<double> x_list;
	vector<pair<double, double>>::iterator it2;
	it = lower_vertex.begin()++;
	it2 = upper_vertex.begin()++;
	x_list1.push_back(temp_min);
	while (it != lower_vertex.end() || it2 != upper_vertex.end()) {
		if (it == lower_vertex.end())
		{
			while (it2 != upper_vertex.end())
			{
				x_list1.push_back(it2->first);
				it2++;
			}
			break;
		}
		else if (it2 == upper_vertex.end()) {
			while (it != lower_vertex.end()) {
				x_list1.push_back(it->first);
				it++;
			}
			break;
		}
		else {
			if (it->first > it2->first) {
				x_list1.push_back(it2->first);
				it2++;
			}
			else if (it->first < it2->first) {
				x_list1.push_back(it->first);
				it++;
			}
			else {
				x_list1.push_back(it->first);
				it++;
				it2++;
			}

		}

	}
	vector<double>::iterator it3;
	vector<double>::iterator it4;
	for (it3 = x_list1.begin(); it3 != x_list1.end(); it3++) {
		if (*it3 - width>temp_min)
			x_list2.push_back(*it3 - width);
	}
	if (x_list2.size() != 0) {
		it3 = x_list1.begin();
		it4 = x_list2.begin();

		if (*it3 < *it4)
		{
			x_list.push_back(*it3);
			it3++;
		}
		else if (*it4 < *it3)
		{
			x_list.push_back(*it4);
			it4++;
		}
		else {
			x_list.push_back(*it3);
			it3++;
			it4++;
		}

		while (it3 != x_list1.end() || it4 != x_list2.end()) {
			if (it3 == x_list1.end())
			{
				while (it4 != x_list2.end())
				{
					if (x_list[x_list.size() - 1] != *it4)
						x_list.push_back(*it4);
					it4++;
				}
				break;
			}
			else if (it4 == x_list2.end()) {
				while (it3 != x_list1.end()) {
					if (x_list[x_list.size() - 1] != *it3)
						x_list.push_back(*it3);
					it3++;
				}
			}
			else {
				if (*it3 > *it4) {
					if (x_list[x_list.size() - 1] != *it4)
						x_list.push_back(*it4);
					it4++;
				}
				else if (*it3 < *it4) {
					if (x_list[x_list.size() - 1] != *it3)
						x_list.push_back(*it3);
					it3++;
				}
				else {
					if (x_list[x_list.size() - 1] != *it3)
						x_list.push_back(*it3);
					it3++;
					it4++;
				}

			}
		}
	}
	else {
		for (it3 = x_list1.begin(); it3 != x_list1.end(); it3++)
			x_list.push_back(*it3);
	}
	int upper_idx = 0;
	int lower_idx = 0;
	int upper_idx2 = 0;
	int lower_idx2 = 0;
	double temp_area = 0;
	double max_area = 0;
	
	vector<partial_area> area_vec;
	queue<partial_area> area_q;
	for (i = 0; i < x_list.size() - 1; i++) {
		if(x_list[i]!=x_list[i+1])
			area_vec.push_back(partial_area(upper_vertex[upper_idx2], upper_vertex[upper_idx2 + 1], lower_vertex[lower_idx2], lower_vertex[lower_idx2 + 1], x_list[i], x_list[i + 1]));
		if (x_list[i + 1] >= upper_vertex[upper_idx2 + 1].first)
		{
			if (abs(upper_vertex[upper_idx2 + 1].first - upper_vertex[upper_idx2].first)<=DBL_EPSILON)
				upper_idx2++;
			upper_idx2++;

		}
		
		if (x_list[i + 1] >= lower_vertex[lower_idx2 + 1].first)
		{
			if (abs(lower_vertex[lower_idx2 + 1].first - lower_vertex[lower_idx2].first)<=DBL_EPSILON)
				lower_idx2++;
			lower_idx2++;


		}

	}
	upper_idx2 = 0;
	lower_idx2 = 0;
	int right_idx = 0;
	double remain = 0;

	if (lower_vertex[lower_vertex.size() - 1].first - temp_min < width)
	{
		for (int i = 0; i < area_vec.size(); i++)
			temp_area += area_vec[i].area;
		return temp_area;
	}

	while (lower_vertex[lower_idx2 + 1].first - temp_min < width)
	{
		lower_idx2++;
		if (lower_idx2 == lower_vertex.size()-1)
			break;
	}
	while (upper_vertex[upper_idx2 + 1].first - temp_min < width )
	{
		upper_idx2++;
		if (upper_idx2 == upper_vertex.size()-1)
			break;
	}
	while (x_list[right_idx + 1] - temp_min < width)
	{
		right_idx++;
		temp_area += area_vec[right_idx].area;
		if (right_idx == area_vec.size() - 1)
		{
			break;
		}
	}
	temp_area -= area_vec[right_idx].area;
	remain = width - (x_list[right_idx] - x_list[1]);
	for (i = 0; i < x_list.size() - 1; i++) {
		if (right_idx == area_vec.size())
			break;
		remain = width - (x_list[right_idx] - x_list[i + 1]);
		double test = area_vec[i].find_max_area(area_vec[right_idx],remain);
		max_area = max(max_area, temp_area + area_vec[i].find_max_area(area_vec[right_idx], remain));
		temp_area -= area_vec[i + 1].area;
		if (x_list[i + 1] + width >= x_list[right_idx + 1])
		{
			temp_area += area_vec[right_idx].area;
			right_idx++;
		}
		
	}
	return max_area;
}

int binary_search(vector<double> y_list, vertex* top, vertex* bot, int start, int end, double width, double height) {
	if (start == end)
		return start;
	int mid1 = (start + end) / 2;
	int mid2 = (start + end) / 2 + 1;
	double y1 = y_list[mid1];
	double y2 = y_list[mid2];
	double u_line= find_max_area(y2, top, bot, width, height);
	double l_line = find_max_area(y1, top, bot, width, height);
	
	if (start + 1 == end)
		if (u_line >= l_line)
			return mid2;
		else
			return mid1;
	else if (u_line >= l_line)
		return binary_search(y_list, top, bot, (start + end) / 2 + 1, end, width, height);
	else if (u_line < l_line)
		return binary_search(y_list, top, bot, (start), (start + end) / 2, width, height);

}

void trans_polygon(vertex* left,vertex* right, vertex** t_top, vertex** t_bot) {
	vertex* temp = left;
	vertex* temp2;
	edge* e_temp;
	*t_bot = new vertex(left->v_y, left->v_x);
	(*t_bot)->edge_ccw = new edge();
	(*t_bot)->edge_ccw->from = *t_bot;
	e_temp = (*t_bot)->edge_ccw;
	e_temp->from = (*t_bot);
	temp = temp->edge_cw->from;
	while (temp != left) {
		e_temp->to = new vertex(temp->v_y, temp->v_x);
		e_temp->to->edge_cw = e_temp;
		if (temp == right)
			*t_top = e_temp->to;
		temp2 = e_temp->to;
		temp2->edge_ccw = new edge();
		temp2->edge_ccw->from = temp2;
		e_temp = temp2->edge_ccw;
		temp = temp->edge_cw->from;
	}
	e_temp->to = *t_bot;
	(*t_bot)->edge_cw = e_temp;

	
}

void find_max_area_side(double x_2, double xy, double y_2, double x_1, double y_1, double constant, vector<pair<double, double>> domain) {
	pair<double, double> pos_temp;

	double temp_max = 0;
	double x_start, x_end,y_start, y_end;
	double a, b;
	double x_max_point;
	double y_max_point;
	for (int i = 0; i < domain.size(); i++) {
		if (i != domain.size() - 1) {
			x_start = domain[i].first;
			x_end = domain[i + 1].first;
			y_start = domain[i].second;
			y_end = domain[i + 1].second;
		}
		else {
			x_start = domain[i].first;
			x_end = domain[0].first;
			y_start = domain[i].second;
			y_end = domain[0].second;
		}

		if (abs(x_start - x_end)<=DBL_EPSILON) {
			y_max_point = -(xy*x_start + y_1) / (2 * y_2);
			x_max_point = x_start;

			temp_max = x_2 * x_start*x_start + xy * x_start*y_max_point + y_2 * y_max_point*y_max_point + x_1 * x_start + y_1 * y_max_point + constant;

			if (global_max < temp_max && y_max_point > min(y_start,y_end) && y_max_point <max(y_start, y_end) ) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}
			y_max_point = y_end;
			temp_max = x_2 * x_start*x_start + xy * x_start*y_max_point + y_2 * y_max_point*y_max_point + x_1 * x_start + y_1 * y_max_point + constant;
			
			if (global_max < temp_max) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}

			y_max_point = y_start;
			temp_max = x_2 * x_start*x_start + xy * x_start*y_max_point + y_2 * y_max_point*y_max_point + x_1 * x_start + y_1 * y_max_point + constant;

			if (global_max < temp_max) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}

		}
		else{
			if (x_start > x_end) {
				swap(x_start, x_end);
				swap(y_start, y_end);
			}
			a = (y_end - y_start) / (x_end - x_start);
			b = y_end - a * x_end;
			
			x_max_point = -(xy*b + 2 * a*b*y_2 + x_1 + y_1 * a) / (2 * x_2 + a * xy + x_2 * xy + 2 * y_2*a*a);
			y_max_point = a * x_max_point + b;

			temp_max = x_2 * x_max_point*x_max_point + xy * x_max_point*(a*x_max_point + b) + y_2 * (a*x_max_point + b)*(a*x_max_point + b) + x_1 * x_max_point + y_1 * (a*x_max_point + b) + constant;
			if (global_max < temp_max && x_max_point < x_end && x_max_point>x_max_point) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}

			x_max_point = x_start;
			y_max_point = a * x_max_point + b;
			temp_max = x_2 * x_start*x_start + xy * x_start*(a*x_start + b) + y_2 * (a*x_start + b)*(a*x_start + b) + x_1 * x_start + y_1 * (a*x_start + b) + constant;
			if (global_max < temp_max) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}

			x_max_point = x_end;
			y_max_point = a * x_max_point + b;
			temp_max = x_2 * x_end*x_end + xy * x_end*(a*x_end + b) + y_2 * (a*x_end + b)*(a*x_end + b) + x_1 * x_end + y_1 * (a*x_end + b) + constant;
			if (global_max < temp_max) {
				global_max = temp_max;
				max_pos.first = x_max_point;
				max_pos.second = y_max_point;
			}
			
		}
		
	}
	
}