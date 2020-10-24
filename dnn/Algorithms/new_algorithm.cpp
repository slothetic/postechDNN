#include<iostream>
#include<algorithm>
#include<vector>
#include<cfloat>
#include<ctime>
using namespace std;
class coeff_set {
public:
	int quad;
	double x2;
	double x1;
	double xy;
	double y2;
	double y1;
	double constant;
	double x_peak;
	double y_peak;
	coeff_set() {
		x2 = x1 = xy = y2 = y1 = constant = 0.0;
	}
	double cal_area(double i_x,double i_y) {
		if ((quad == 1 || quad == 4) && i_x > x_peak)
			i_x = x_peak;
		else if ((quad == 2 || quad == 3) && i_x < x_peak)
			i_x = x_peak;
		if ((quad == 1 || quad == 2) && i_y > y_peak)
			i_y = y_peak;
		else if ((quad == 3 || quad == 4) && i_y < y_peak)
			i_y = y_peak;
		return x2 * i_x*i_x + y2 * i_y*i_y + xy * i_x*i_y + x1 * i_x + y1 * i_y + constant;
	}
};

class vertex {
public:
	class edge* edge_cw;
	class edge* edge_ccw;
	double v_x;
	double v_y;
	vertex(double in_x, double in_y) { v_x = in_x; v_y = in_y; }
};

class edge {
public:
	vertex * from;
	vertex* to;
	double temp_x1;
	double temp_x2;
	double temp_y1;
	double temp_y2;
	pair<bool, double> gradient; // if gradient.first equals to true, edge is perpendicular to x-axis
	
	double x_intersect(double line_y) {
		if (abs(from->v_x - to->v_x) <= DBL_EPSILON)
			return from->v_x;
		else {
			return from->v_x + (from->v_y - line_y)*(to->v_x - from->v_x) / (from->v_y - to->v_y);
		}
	}
	double y_intersect(double line_x) {
		if (abs(from->v_y - to->v_y) <= DBL_EPSILON)
			return from->v_y;
		else
			return from->v_y + (line_x - from->v_x)*(from->v_y - to->v_y) / (from->v_x - to->v_x);
	}
	pair<bool, pair<double, double>> y_check(double line_y) {
		double y_high = max(from->v_y, to->v_y);
		double y_low = min(from->v_y, to->v_y);
		if (y_high > line_y && line_y >= y_low && (y_high != y_low)) {
			return pair<bool, pair<double, double>>(true, pair<double, double>(x_intersect(line_y), line_y));
		}
		else {
			return pair<bool, pair<double, double>>(false, pair<double, double>(0, 0));
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
			return pair<bool, pair<double, double>>(true, pair<double, double>(line_x, y_intersect(line_x)));
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

class area {
public:
	pair<double, double> from, to;
	double max_area;
	double accumulated_area;
	double coeff_2;
	double coeff_1;
	double coeff_con;
	double partial_area(double k) {
		return coeff_2 * k*k + coeff_1 * k + coeff_con;
	}
	area() {
		max_area = 0;
		coeff_2 = 0;
		coeff_1 = 0;
		coeff_con = 0;
		accumulated_area = 0;
	}
};

class area_struct {
public:
	int quadrant;
	vector<pair<double, double>> vertices;
	double x_peak, y_peak;
	double y_of_x_peak;
	double x_of_y_peak;
	double empty;
	vector<area> vertical_area1;
	vector<area> horizon_area1;
	vector<area> vertical_area2;
	vector<area> horizon_area2;
	void set_area();
	double cal_area(double input_x, double input_y, int v1_idx = -2, int v2_idx = -2, int h1_idx = -2, int h2_idx = -2);
	int find_index(double value, int type);
	void cal_coeff(double input_x, double input_y, double* coeff_x2, double* coeff_x1, double* coeff_y2, double* coeff_y1, double* coeff_xy,double* coeff_con);
};
int x_binary_search(vector<double>& x_list, vector<double>& y_list, int start, int end);
int y_binary_search(vector<double>& x_list, vector<double>& y_list, int start, int end);
double binary_search_horizon(vector<double>& x_list, double y_pos, int start, int end,int** coeff_list);
double binary_search_vertical(double x_pos, vector<double>& y_list, int start, int end,int** coeff_list);
double sum_area(double pos_x, double pos_y,int** coeff);
pair<double, double>find_intersect(pair<double, double> a, pair<double, double> b, double axis, int type);
int binary_search_index(vector<area>& vec, double value, int type, int start, int end, int order);
bool check_external(pair<double, double> from, pair<double, double> to, pair<double, double> point);
double cal_area_coeff(double input_x, double input_y, double c_x2, double c_xy, double c_y2, double c_x1, double c_y1, double c_con);
void find_max_point(double x_low, double x_high, double y_low, double y_high, vertex* top, vertex* bot, vertex* left, vertex* right);
int edge_check(edge* e_temp, double left, double right, double top, double bot);
bool recursive_find(double x_low, double x_high, double y_low, double y_high, double width, double height, vector<int> e_check_list, edge* e_list[], int d_num, vector<pair<double, double>> domain);
//double cal_area_coeff(double input_x, double input_y, double c_x2, double c_xy, double c_y2, double c_x1, double c_y1, double c_con);
pair<double, double> line_intersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
pair<double, double> make_point(vector<pair<double, double>> domain);
void find_max_area_side(double x_2, double xy, double y_2, double x_1, double y_1, double constant, vector<pair<double, double>> domain);
bool check_point(double p_x, double p_y, vector<pair<double, double>>domain);
int check_shift(area_struct& quad, int idx, double value, int type);

vertex* top_vertex;
vertex* bot_vertex;
vertex* left_vertex;
vertex* right_vertex;
double x_axis, y_axis;
area_struct quad1, quad2, quad3, quad4;
double width, height;
double global_max;
pair<double, double> max_pos;

int main() {

	int num_v;
	int phase;
	
	quad1.quadrant = 1;
	quad2.quadrant = 2;
	quad3.quadrant = 3;
	quad4.quadrant = 4;
	
	double left_side, right_side, top_side, bottom_side;
	vertex* v_temp;
	vertex* start_node;
	vertex* t_bot = NULL;
	vertex* t_top = NULL;
	
	cout << "width and heigth:";
	cin >> width >> height;
	double temp_x, temp_y;
	double x_domain_min;
	double y_domain_min;
	double height_max;
	double width_max;
	pair<double, double> temp_vertex;
	cin >> num_v;
	cin >> temp_x >> temp_y;//아무 vertex를 시작으로 반시계 방향 순서로 입력
	left_side = right_side = temp_x;
	top_side = bottom_side = temp_y;

	v_temp = new vertex(temp_x, temp_y);
	start_node = v_temp;
	top_vertex = v_temp;
	bot_vertex = v_temp;
	left_vertex = v_temp;
	right_vertex = v_temp;

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
		if (i == num_v - 1) {
			v_temp->edge_ccw = new edge();
			v_temp->edge_ccw->from = v_temp;
			v_temp->edge_ccw->to = start_node;
			v_temp->edge_ccw->set_gradient();
			start_node->edge_cw = v_temp->edge_ccw;
		}

	}

	clock_t start = clock();
	//////////find position which has maximum height and width
	v_temp = left_vertex;
	pair<bool, pair<double, double>> check_temp;
	while (v_temp != right_vertex) {
		
		check_temp=v_temp->edge_ccw->x_check(top_vertex->v_x);
		if (check_temp.first == true) {
			height_max = top_vertex->v_y - check_temp.second.second;
			y_axis = check_temp.second.first;
			start_node = top_vertex;
			break;
		}
		v_temp = v_temp->edge_ccw->to;
	}

	v_temp = left_vertex;
	while (v_temp != right_vertex) {
		check_temp = v_temp->edge_cw->x_check(bot_vertex->v_x);
		if (check_temp.first == true) {
			if (height_max < check_temp.second.second - bot_vertex->v_y)
			{
				y_axis = check_temp.second.first;
				start_node = bot_vertex;
				break;
			}
		}
		v_temp = v_temp->edge_cw->from;
	}
	
	v_temp = bot_vertex;
	while (v_temp != top_vertex) {
		check_temp = v_temp->edge_ccw->y_check(left_vertex->v_y);
		if (check_temp.first == true) {
			width_max = check_temp.second.first - left_vertex->v_x;
			x_axis = check_temp.second.second;
			break;
		}
		v_temp = v_temp->edge_ccw->to;
	}

	v_temp = bot_vertex;
	while (v_temp != top_vertex) {
		check_temp = v_temp->edge_cw->y_check(right_vertex->v_y);
		if (check_temp.first == true) {
			if (width_max < right_vertex->v_x - check_temp.second.first) {
				x_axis = check_temp.second.second;
				break;
			}
		}
		v_temp = v_temp->edge_cw->from;
	}

	cout << x_axis << "//" << y_axis << endl;
	////////////////////////////////////
	bool vertex_check = 0;
	////push_back vertices to each quadrant
	if (start_node == bot_vertex) {
		
		phase = 4;
		while (!(phase == 3 && start_node == bot_vertex)) {
			
			switch (phase) {
			case 1:
				quad1.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 2:
				quad2.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 3:
				quad3.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 4:
				quad4.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			}

			start_node = start_node->edge_ccw->to;
			if (phase == 4 && start_node->v_y >= x_axis)
			{
				temp_vertex = find_intersect(quad4.vertices[quad4.vertices.size() - 1], pair < double, double >(start_node->v_x, start_node->v_y), x_axis, 2);
				quad4.vertices.push_back(temp_vertex);
				if (start_node->v_y != x_axis)
					quad1.vertices.push_back(temp_vertex);
				phase = 1;
			}
			else if (phase == 1 && start_node->v_x <= y_axis) {
				temp_vertex = find_intersect(quad1.vertices[quad1.vertices.size() - 1], pair < double, double >(start_node->v_x, start_node->v_y), y_axis, 1);
				quad1.vertices.push_back(temp_vertex);
				if (start_node->v_x != y_axis)
					quad2.vertices.push_back(temp_vertex);
				phase = 2;
			}
			else if (phase == 2 && start_node->v_y <= x_axis) {
				temp_vertex = find_intersect(quad2.vertices[quad2.vertices.size() - 1], pair < double, double >(start_node->v_x, start_node->v_y), x_axis, 2);
				quad2.vertices.push_back(temp_vertex);
				if (start_node->v_y != x_axis)
					quad3.vertices.push_back(temp_vertex);
				phase = 3;
			}
		}
		quad3.vertices.push_back(pair<double, double>(bot_vertex->v_x, bot_vertex->v_y));

	}
	else {
		phase = 2;
		while (!(phase == 1 && start_node == top_vertex)) {
			switch (phase) {
			case 1:
				quad1.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 2:
				quad2.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 3:
				quad3.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			case 4:
				quad4.vertices.push_back(pair<double, double>(start_node->v_x, start_node->v_y));
				break;
			}
			start_node = start_node->edge_ccw->to;
			if (phase == 2 && start_node->v_y <= x_axis) {
				temp_vertex = find_intersect(quad2.vertices[quad2.vertices.size() - 1], pair < double,double >(start_node->v_x,start_node->v_y) , x_axis, 2);
				quad2.vertices.push_back(temp_vertex);
				if (start_node->v_y != x_axis)
					quad3.vertices.push_back(temp_vertex);
				phase = 3;
			}
			else if (phase == 3 && start_node->v_x >= y_axis) {
				temp_vertex = find_intersect(quad3.vertices[quad3.vertices.size() - 1], pair < double, double >(start_node->v_x, start_node->v_y), y_axis, 1);
				quad3.vertices.push_back(temp_vertex);
				if (start_node->v_x != y_axis)
					quad4.vertices.push_back(temp_vertex);
				phase = 4;
			}
			else if (phase == 4 && start_node->v_y >= x_axis) {
				temp_vertex = find_intersect(quad4.vertices[quad4.vertices.size() - 1], pair < double, double >(start_node->v_x, start_node->v_y), x_axis, 2);
				quad4.vertices.push_back(temp_vertex);
				if (start_node->v_y != x_axis)
					quad1.vertices.push_back(temp_vertex);
				phase = 1;
			}
		}
		quad1.vertices.push_back(pair<double, double>(top_vertex->v_x, top_vertex->v_y));
	}


	/// set the value x_peak and y_peak for each area_struct
	if (right_vertex->v_y > x_axis) {
		quad1.x_peak = right_vertex->v_x;
		quad1.y_of_x_peak = right_vertex->v_y;
		quad4.x_peak = quad4.vertices[quad4.vertices.size() - 1].first;
		quad4.y_of_x_peak = quad4.vertices[quad4.vertices.size() - 1].second;
	}
	else if (right_vertex->v_y < x_axis) {
		quad1.x_peak = quad1.vertices[0].first;
		quad1.y_of_x_peak = quad1.vertices[0].second;
		quad4.x_peak = right_vertex->v_x;
		quad4.y_of_x_peak = right_vertex->v_y;
	}
	else {
		quad1.x_peak = right_vertex->v_x;
		quad1.y_of_x_peak = right_vertex->v_y;
		quad4.x_peak = right_vertex->v_x;
		quad4.y_of_x_peak = right_vertex->v_y;
	}

	if (left_vertex->v_y > x_axis) {
		quad2.x_peak = left_vertex->v_x;
		quad2.y_of_x_peak = left_vertex->v_y;
		quad3.x_peak = quad3.vertices[0].first;
		quad3.y_of_x_peak = quad3.vertices[0].second;
	}
	else if (left_vertex->v_y < x_axis) {
		quad2.x_peak = quad2.vertices[quad2.vertices.size() - 1].first;
		quad2.y_of_x_peak = quad2.vertices[quad2.vertices.size() - 1].second;
		quad3.x_peak = left_vertex->v_x;
		quad3.y_of_x_peak = left_vertex->v_y;
	}
	else {
		quad2.x_peak = left_vertex->v_x;
		quad2.y_of_x_peak = left_vertex->v_y;
		quad3.x_peak = left_vertex->v_x;
		quad3.y_of_x_peak = left_vertex->v_y;
	}

	if (top_vertex->v_x < y_axis) {
		quad2.y_peak = top_vertex->v_y;
		quad2.x_of_y_peak = top_vertex->v_x;
		quad1.y_peak = quad1.vertices[quad1.vertices.size() - 1].second;
		quad2.x_of_y_peak = quad1.vertices[quad1.vertices.size() - 1].first;
	}
	else if (top_vertex->v_x > y_axis) {
		quad1.y_peak = top_vertex->v_y;
		quad1.x_of_y_peak = top_vertex->v_x;
		quad2.y_peak = quad2.vertices[0].second;
		quad2.x_of_y_peak = quad2.vertices[0].first;
	}
	else {
		quad1.y_peak = top_vertex->v_y;
		quad1.x_of_y_peak = top_vertex->v_x;
		quad2.y_peak = top_vertex->v_y;
		quad2.x_of_y_peak = top_vertex->v_x;
	}
	
	if (bot_vertex->v_x < y_axis) {
		quad3.y_peak = bot_vertex->v_y;
		quad3.x_of_y_peak = bot_vertex->v_x;
		quad4.y_peak = quad4.vertices[0].second;
		quad4.x_of_y_peak = quad4.vertices[0].first;
	}
	else if (bot_vertex->v_x > y_axis) {
		quad3.y_peak = quad3.vertices[quad3.vertices.size() - 1].second;
		quad3.x_of_y_peak = quad3.vertices[quad3.vertices.size() - 1].first;
		quad4.y_peak = bot_vertex->v_y;
		quad4.x_of_y_peak = bot_vertex->v_x;
	}
	else {
		quad3.y_peak = bot_vertex->v_y;
		quad3.x_of_y_peak = bot_vertex->v_x;
		quad4.y_peak = bot_vertex->v_y;
		quad4.x_of_y_peak = bot_vertex->v_x;
	}
	
	
	//set up the area functions for each quadrant
	quad1.set_area();
	quad2.set_area();
	quad3.set_area();
	quad4.set_area();
	///////////


	/* test code
	double t_x, t_y;
	double t_w, t_h;
	t_x = 0.5;
	t_y = 1.81047;
	t_w = 16.5;
	t_h = 3.09;
	cout << quad1.cal_area(t_x + t_w, t_y + t_h) <<" + "<< quad2.cal_area(t_x, t_y + t_h) << " + " << quad3.cal_area(t_x, t_y) << " + " << quad4.cal_area(t_x + t_w, t_y)<<endl;
	cout << quad1.cal_area(t_x+t_w, t_y+t_h) + quad2.cal_area(t_x, t_y+t_h) + quad3.cal_area(t_x, t_y) + quad4.cal_area(t_x+t_w, t_y);
	*/

	/// create x_list and y_list
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
	vector<double>::iterator it, it2;
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
	//////

	////minimize the domain of x and y. 
	x_list.clear();
	it = x_list_res.begin();
	x_domain_min = max(y_axis - width, left_side);
	x_list.push_back(x_domain_min);
	while (it != x_list_res.end()) {
		if (x_domain_min < *it && abs(x_domain_min - *it) >= DBL_EPSILON) {
			if (*it < y_axis && abs(*it - y_axis) >= DBL_EPSILON) {
				x_list.push_back(*it);
			}
		}
		it++;
	}
	x_list.push_back(y_axis);

	y_list.clear();
	it = y_list_res.begin();
	y_domain_min = max(x_axis - height, bottom_side);
	y_list.push_back(y_domain_min);
	while (it != y_list_res.end()) {
		if (y_domain_min < *it && abs(y_domain_min - *it) >= DBL_EPSILON) {
			if (*it < x_axis && abs(*it - x_axis) >= DBL_EPSILON) {
				y_list.push_back(*it);
			}
		}
		it++;
	}
	y_list.push_back(x_axis);
	/////
	int x_idx, y_idx;
	
	x_idx = x_binary_search(x_list, y_list, 0, x_list.size() - 1);
	y_idx = y_binary_search(x_list, y_list, 0, y_list.size() - 1);
	cout << "x: " << x_list[x_idx] << endl;
	cout << "y: " << y_list[y_idx] << endl;
	
	if (x_idx != 0) {
		if (y_idx != 0)
			find_max_point(x_list[x_idx - 1], x_list[x_idx], y_list[y_idx - 1], y_list[y_idx], top_vertex, bot_vertex, left_vertex, right_vertex);
		if (y_idx != y_list.size() - 1)
			find_max_point(x_list[x_idx - 1], x_list[x_idx], y_list[y_idx], y_list[y_idx + 1], top_vertex, bot_vertex, left_vertex, right_vertex);
	}

	if (x_idx != x_list.size() - 1 && y_idx != 0)
		find_max_point(x_list[x_idx], x_list[x_idx + 1], y_list[y_idx - 1], y_list[y_idx], top_vertex, bot_vertex, left_vertex, right_vertex);
	if (x_idx != x_list.size() - 1 && y_idx != y_list.size() - 1)
		find_max_point(x_list[x_idx], x_list[x_idx + 1], y_list[y_idx], y_list[y_idx + 1], top_vertex, bot_vertex, left_vertex, right_vertex);
	cout.precision(8);
	cout << "answer: (x,y) = (" << max_pos.first << " , " << max_pos.second << " ) /// area: " << global_max << endl;
	cout.precision(8);
	clock_t end = clock();
	cout << "\nTime: " << (double)(end - start) / CLOCKS_PER_SEC << endl;
	return 0;
}

int check_shift(area_struct& quad, int idx, double value,int type) {
	
	if (idx == -1) {
		return quad.find_index(value, type);
	}
	
	if (type == 1) {
		if (max(quad.vertical_area1[idx].from.first, quad.vertical_area1[idx].to.first) + DBL_EPSILON >= value && value >= min(quad.vertical_area1[idx].from.first, quad.vertical_area1[idx].to.first) + DBL_EPSILON) {
			
			return idx;
		}
	}
	else if(type==2){
		if (max(quad.vertical_area2[idx].from.first, quad.vertical_area2[idx].to.first) + DBL_EPSILON >= value && value >= min(quad.vertical_area2[idx].from.first, quad.vertical_area2[idx].to.first) + DBL_EPSILON) {
			
			return idx;
		}
	}
	else if (type == 3) {
		if (max(quad.horizon_area1[idx].from.second, quad.horizon_area1[idx].to.second) + DBL_EPSILON >= value && value >= min(quad.horizon_area1[idx].from.second, quad.horizon_area1[idx].to.second) + DBL_EPSILON) {
			
			return idx;
		}
	}
	else if (type == 4) {
		if (max(quad.horizon_area2[idx].from.second, quad.horizon_area2[idx].to.second) + DBL_EPSILON >= value && value >= min(quad.horizon_area2[idx].from.second, quad.horizon_area2[idx].to.second) + DBL_EPSILON) {
			
			return idx;
		}
	}
	return quad.find_index(value, type);
}

int x_binary_search(vector<double>& x_list, vector<double>& y_list,int start, int end) {
	if (start == end)
		return start;
	int mid = (start + end) / 2;
	int ** coeff_list = new int*[4];
	for (int i = 0; i < 4; i++) {
		coeff_list[i] = new int[4];
		for (int j = 0; j < 4; j++)
			coeff_list[i][j] = -2;
	}
	coeff_list[0][0]=quad1.find_index(x_list[mid]+width, 1);
	coeff_list[0][1] = quad1.find_index(x_list[mid]+width, 2);
	cout << coeff_list[0][1] << endl;
	coeff_list[1][0] = quad2.find_index(x_list[mid], 1);
	coeff_list[1][1] = quad2.find_index(x_list[mid], 2);
	coeff_list[2][0] = quad3.find_index(x_list[mid] , 1);
	coeff_list[2][1] = quad3.find_index(x_list[mid] , 2);
	coeff_list[3][0] = quad4.find_index(x_list[mid]+width, 1);
	coeff_list[3][1] = quad4.find_index(x_list[mid]+width, 2);
	double area1 = binary_search_vertical(x_list[mid], y_list, 0, y_list.size() - 1,coeff_list);
	
	coeff_list[0][0] = check_shift(quad1, coeff_list[0][0], x_list[mid + 1] + width, 1);
	coeff_list[0][1] = check_shift(quad1,coeff_list[0][1],x_list[mid + 1] + width, 2);
	coeff_list[1][0] = check_shift(quad2,coeff_list[1][0],x_list[mid + 1], 1);
	coeff_list[1][1] = check_shift(quad2,coeff_list[1][1],x_list[mid + 1], 2);
	coeff_list[2][0] = check_shift(quad3,coeff_list[2][0],x_list[mid + 1], 1);
	coeff_list[2][1] = check_shift(quad3,coeff_list[2][1],x_list[mid + 1], 2);
	coeff_list[3][0] = check_shift(quad4,coeff_list[3][0],x_list[mid + 1] + width, 1);
	coeff_list[3][1] = check_shift(quad4, coeff_list[3][1], x_list[mid + 1] + width, 2);
	/*
	coeff_list[0][0] = quad1.find_index(x_list[mid+1] + width, 1);
	coeff_list[0][1] = quad1.find_index(x_list[mid+1] + width, 2);
	coeff_list[1][0] = quad2.find_index(x_list[mid+1], 1);
	coeff_list[1][1] = quad2.find_index(x_list[mid+1], 2);
	coeff_list[2][0] = quad3.find_index(x_list[mid+1], 1);
	coeff_list[2][1] = quad3.find_index(x_list[mid+1], 2);
	coeff_list[3][0] = quad4.find_index(x_list[mid+1] + width, 1);
	coeff_list[3][1] = quad4.find_index(x_list[mid+1] + width, 2);
	*/
	double area2 = binary_search_vertical(x_list[mid + 1], y_list, 0, y_list.size() - 1,coeff_list);
	delete[] coeff_list[0];
	delete[] coeff_list[1];
	delete[] coeff_list[2];
	delete[] coeff_list[3];
	delete[] coeff_list;
	
	if (area1 >= area2)
		return x_binary_search(x_list, y_list, start, mid);
	else
		return x_binary_search(x_list, y_list, mid + 1, end);
}

int y_binary_search(vector<double>& x_list, vector<double>& y_list, int start, int end) {
	if (start == end)
		return start;
	int mid = (start + end) / 2;
	int ** coeff_list = new int*[4];
	for (int i = 0; i < 4; i++) {
		coeff_list[i] = new int[4];
		for (int j = 0; j < 4; j++)
			coeff_list[i][j] = -2;
	}
	coeff_list[0][2] = quad1.find_index(y_list[mid] + height, 3);
	coeff_list[0][3] = quad1.find_index(y_list[mid] + height, 4);
	coeff_list[1][2] = quad2.find_index(y_list[mid]+height, 3);
	coeff_list[1][3] = quad2.find_index(y_list[mid]+height, 4);
	coeff_list[2][2] = quad3.find_index(y_list[mid], 3);
	coeff_list[2][3] = quad3.find_index(y_list[mid], 4);
	coeff_list[3][2] = quad4.find_index(y_list[mid] , 3);
	coeff_list[3][3] = quad4.find_index(y_list[mid] , 4);
	double area1 = binary_search_horizon(x_list, y_list[mid], 0, x_list.size() - 1,coeff_list);
	coeff_list[0][2] = check_shift(quad1,coeff_list[0][2],y_list[mid + 1] + height, 3);
	coeff_list[0][3] = check_shift(quad1,coeff_list[0][3],y_list[mid + 1] + height, 4);
	coeff_list[1][2] = check_shift(quad2,coeff_list[1][2],y_list[mid + 1]+height, 3);
	coeff_list[1][3] = check_shift(quad2,coeff_list[1][3],y_list[mid + 1]+height, 4);
	coeff_list[2][2] = check_shift(quad3,coeff_list[2][2],y_list[mid + 1], 3);
	coeff_list[2][3] = check_shift(quad3,coeff_list[2][3],y_list[mid + 1], 4);
	coeff_list[3][2] = check_shift(quad4,coeff_list[3][2],y_list[mid + 1] , 3);
	coeff_list[3][3] = check_shift(quad4,coeff_list[3][3],y_list[mid + 1] , 4);
	
	double area2 = binary_search_horizon(x_list, y_list[mid + 1], 0, x_list.size() - 1,coeff_list);
	delete[] coeff_list[0];
	delete[] coeff_list[1];
	delete[] coeff_list[2];
	delete[] coeff_list[3];
	delete[] coeff_list;
	if (area1 >= area2)
		return y_binary_search(x_list, y_list, start, mid);
	else
		return y_binary_search(x_list, y_list, mid + 1, end);
}

double binary_search_horizon(vector<double>& x_list,double y_pos, int start, int end,int** coeff_list) {
	if (start == end) {
		vector<double> x_list_partial;
		double coeff_x2 = 0.0;
		double coeff_x1 = 0.0;
		double coeff_y2 = 0.0;
		double coeff_y1 = 0.0;
		double coeff_xy = 0.0;
		double coeff_con = 0.0;
		double d_max, d_min;
		double max_area;
		double x_pos;
		d_max = d_min = start;
		int idx_temp;
		pair<double, double> pos_temp;
		if (start != 0) {
			x_list_partial.push_back(x_list[start - 1]);
			d_min = x_list[start - 1];
		}
		x_list_partial.push_back(x_list[start]);
		if (start != x_list.size() - 1) {
			x_list_partial.push_back(x_list[start + 1]);
			d_max = x_list[start + 1];
		}

		//find intersection
		idx_temp = coeff_list[0][2];// quad1.find_index(y_pos + height, 3);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad1.horizon_area1[idx_temp].from, quad1.horizon_area1[idx_temp].to, y_pos+height, 2);
			pos_temp.first -= width;
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[0][3];// quad1.find_index(y_pos + height, 4);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad1.horizon_area2[idx_temp].from, quad1.horizon_area2[idx_temp].to, y_pos+height, 2);
			pos_temp.first -= width;
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[1][2];// quad2.find_index(y_pos + height, 3);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad2.horizon_area1[idx_temp].from, quad2.horizon_area1[idx_temp].to, y_pos + height, 2);
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[1][3];// quad2.find_index(y_pos + height, 4);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad2.horizon_area2[idx_temp].from, quad2.horizon_area2[idx_temp].to, y_pos + height, 2);
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[2][2];// quad3.find_index(y_pos, 3);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad3.horizon_area1[idx_temp].from, quad3.horizon_area1[idx_temp].to, y_pos , 2);
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[2][3];// quad3.find_index(y_pos, 4);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad3.horizon_area2[idx_temp].from, quad3.horizon_area2[idx_temp].to, y_pos , 2);
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[3][2];// quad4.find_index(y_pos, 3);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad4.horizon_area1[idx_temp].from, quad4.horizon_area1[idx_temp].to, y_pos, 2);
			pos_temp.first -= width;
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		idx_temp = coeff_list[3][3];// quad4.find_index(y_pos, 4);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad4.horizon_area2[idx_temp].from, quad4.horizon_area2[idx_temp].to, y_pos, 2);
			pos_temp.first -= width;
			if (pos_temp.first < d_max && pos_temp.first>d_min && abs(pos_temp.first - d_max) > DBL_EPSILON && abs(pos_temp.first - d_min) > DBL_EPSILON) {
				x_list_partial.push_back(pos_temp.first);
			}
		}
		/////////////
		sort(x_list_partial.begin(), x_list_partial.end());
		max_area = sum_area(x_list_partial[0],y_pos,coeff_list);
		for (int i = 0; i < x_list_partial.size() - 1; i++) {
			coeff_x2 = 0.0;
			coeff_x1 = 0.0;
			coeff_y2 = 0.0;
			coeff_y1 = 0.0;
			coeff_xy = 0.0;
			coeff_con = 0.0;
			quad1.cal_coeff((x_list_partial[i] + x_list_partial[i + 1]) / 2.0+width,y_pos+ height, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad2.cal_coeff((x_list_partial[i] + x_list_partial[i + 1]) / 2.0, y_pos + height, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad3.cal_coeff((x_list_partial[i] + x_list_partial[i + 1]) / 2.0, y_pos, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad4.cal_coeff((x_list_partial[i] + x_list_partial[i + 1]) / 2.0 + width, y_pos, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			if (abs(coeff_x2) > DBL_EPSILON) {
				x_pos = -(coeff_xy*y_pos + coeff_x1) / (2.0*coeff_x2);
				if (x_pos<d_max && x_pos>d_min || abs(x_pos - d_min) <= DBL_EPSILON || abs(x_pos - d_max) <= DBL_EPSILON) {
					max_area = max(max_area, cal_area_coeff(x_pos, y_pos, coeff_x2, coeff_xy, coeff_y2, coeff_x1, coeff_y1, coeff_con));
					
				}
			}
			max_area = max(max_area, cal_area_coeff(x_list_partial[i + 1], y_pos, coeff_x2, coeff_xy, coeff_y2, coeff_x1, coeff_y1, coeff_con));
			
		}

		return max_area;
	}
	
	int mid = (start + end) / 2;
	double area1 = sum_area(x_list[mid], y_pos,coeff_list);
	double area2 = sum_area(x_list[mid + 1], y_pos,coeff_list);
	if (area1 >= area2) 
		return binary_search_horizon(x_list, y_pos, start, mid,coeff_list);
	else
		return binary_search_horizon(x_list, y_pos, mid + 1, end,coeff_list);
	
}

double binary_search_vertical(double x_pos, vector<double>& y_list, int start, int end,int** coeff_list) {
	if (start == end) {
		vector<double> y_list_partial;
		double coeff_x2 = 0.0;
		double coeff_x1 = 0.0;
		double coeff_y2 = 0.0;
		double coeff_y1 = 0.0;
		double coeff_xy = 0.0;
		double coeff_con = 0.0;
		double d_max, d_min;
		double max_area;
		double temp_area;
		double y_pos;
		d_max = d_min = start;
		int idx_temp;
		pair<double, double> pos_temp;
		if (start != 0) {
			y_list_partial.push_back(y_list[start-1]);
			d_min = y_list[start - 1];
		}
		y_list_partial.push_back(y_list[start]);
		if (start != y_list.size() - 1) {
			y_list_partial.push_back(y_list[start + 1]);
			d_max = y_list[start + 1];
		}

		//find intersection
		idx_temp = quad1.find_index(x_pos + width, 1);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad1.vertical_area1[idx_temp].from, quad1.vertical_area1[idx_temp].to, x_pos + width, 1);
			pos_temp.second -= height;
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad1.find_index(x_pos + width, 2);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad1.vertical_area2[idx_temp].from, quad1.vertical_area2[idx_temp].to, x_pos + width, 1);
			pos_temp.second -= height;
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad2.find_index(x_pos, 1);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad2.vertical_area1[idx_temp].from, quad2.vertical_area1[idx_temp].to, x_pos, 1);
			pos_temp.second -= height;
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad2.find_index(x_pos, 2);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad2.vertical_area2[idx_temp].from, quad2.vertical_area2[idx_temp].to, x_pos, 1);
			pos_temp.second -= height;
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad3.find_index(x_pos, 1);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad3.vertical_area1[idx_temp].from, quad3.vertical_area1[idx_temp].to, x_pos, 1);
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad3.find_index(x_pos, 2);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad3.vertical_area2[idx_temp].from, quad3.vertical_area2[idx_temp].to, x_pos, 1);
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad4.find_index(x_pos + width, 1);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad4.vertical_area1[idx_temp].from, quad4.vertical_area1[idx_temp].to, x_pos + width, 1);
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}
		idx_temp = quad4.find_index(x_pos + width, 2);
		if (idx_temp != -1) {
			pos_temp = find_intersect(quad4.vertical_area2[idx_temp].from, quad4.vertical_area2[idx_temp].to, x_pos + width, 1);
			if (pos_temp.second < d_max && pos_temp.second>d_min && abs(pos_temp.second - d_max) > DBL_EPSILON && abs(pos_temp.second - d_min) > DBL_EPSILON) {
				y_list_partial.push_back(pos_temp.second);
			}
		}

		/////////////
		sort(y_list_partial.begin(), y_list_partial.end());
		max_area = sum_area(x_pos,y_list_partial[0],coeff_list);
		for (int i = 0; i < y_list_partial.size() - 1; i++) {
			coeff_x2 = 0.0;
			coeff_x1 = 0.0;
			coeff_y2 = 0.0;
			coeff_y1 = 0.0;
			coeff_xy = 0.0;
			coeff_con = 0.0;
			quad1.cal_coeff(x_pos+width, (y_list_partial[i]+y_list_partial[i+1])/2.0+height, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad2.cal_coeff(x_pos, (y_list_partial[i] + y_list_partial[i + 1]) / 2.0 + height, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad3.cal_coeff(x_pos, (y_list_partial[i] + y_list_partial[i + 1]) / 2.0, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			quad4.cal_coeff(x_pos+width, (y_list_partial[i] + y_list_partial[i + 1]) / 2.0, &coeff_x2, &coeff_x1, &coeff_y2, &coeff_y1, &coeff_xy,&coeff_con);
			if (abs(coeff_y2) > DBL_EPSILON) {
				y_pos = -(coeff_xy*x_pos + coeff_y1) / (2.0*coeff_y2);
				if (y_pos<d_max && y_pos>d_min || abs(y_pos - d_min) <= DBL_EPSILON || abs(y_pos - d_max) <= DBL_EPSILON) {
					max_area = max(max_area,cal_area_coeff(x_pos, y_pos,coeff_x2,coeff_xy,coeff_y2,coeff_x1,coeff_y1,coeff_con));
				}
			}
			max_area = max(max_area, cal_area_coeff(x_pos, y_list_partial[i + 1], coeff_x2, coeff_xy, coeff_y2, coeff_x1, coeff_y1, coeff_con));
		}
		
		return max_area;
	}
		
	int mid = (start + end) / 2;
	
	double area1 = sum_area(x_pos, y_list[mid],coeff_list);
	double area2 = sum_area(x_pos, y_list[mid + 1],coeff_list);
	if (area1 >= area2)
		return binary_search_vertical(x_pos, y_list, start, mid, coeff_list);
	else
		return binary_search_vertical(x_pos, y_list, mid + 1, end,coeff_list);
}

double sum_area(double pos_x, double pos_y,int** coeff) {
	return quad1.cal_area(pos_x + width, pos_y + height,coeff[0][0],coeff[0][1],coeff[0][2],coeff[0][3]) + quad2.cal_area(pos_x, pos_y + height,coeff[1][0], coeff[1][1], coeff[1][2], coeff[1][3]) + quad3.cal_area(pos_x, pos_y,coeff[2][0], coeff[2][1], coeff[2][2], coeff[2][3]) + quad4.cal_area(pos_x + width, pos_y, coeff[3][0], coeff[3][1], coeff[3][2], coeff[3][3]);
}

pair<double, double> find_intersect(pair<double, double> a, pair<double, double> b,double axis ,int type) {// type 1: line is vertical / type 2: line is horizontal
	double gradient; 
	if (type == 1) {
		gradient = (b.second - a.second) / (b.first - a.first);
		return pair<double, double>(axis, gradient*(axis - b.first) + b.second);
	}
	else {
		
		if (abs(b.first - a.first) <= FLT_EPSILON)
			return pair<double, double>(a.first, axis);
		else {
			gradient = (b.second - a.second) / (b.first - a.first);
			return pair<double, double>((axis - b.second) / gradient + b.first, axis);
		}
	}
}

void area_struct::set_area() {
	int phase = 1;
	int index;
	double gradient;
	double constant;
	if (quadrant == 1) {
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].first - x_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i+1].second-vertices[i].second)<=DBL_EPSILON) {
				
				if (phase == 1) {
					vertical_area1.push_back(area());
					index = vertical_area1.size() - 1;
					vertical_area1[index].from = vertices[i];
					vertical_area1[index].to = vertices[i + 1];
				}
				else {
					vertical_area2.push_back(area());
					index = vertical_area2.size() - 1;
					vertical_area2[index].to = vertices[i + 1];
					vertical_area2[index].from = vertices[i];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				vertical_area1.push_back(area());
				index = vertical_area1.size() - 1;
				vertical_area1[index].max_area = (vertices[i].second + vertices[i + 1].second - 2 * x_axis)*(vertices[i + 1].first - vertices[i].first) / 2.0;
				vertical_area1[index].from = vertices[i];
				vertical_area1[index].to = vertices[i + 1];
				vertical_area1[index].coeff_2 = -gradient / 2.0;
				vertical_area1[index].coeff_1 = (gradient * vertices[i + 1].first-constant -  vertices[i + 1].second + 2 * x_axis) / 2.0;
				vertical_area1[index].coeff_con = (constant*vertices[i + 1].first + vertices[i + 1].first*vertices[i + 1].second - 2* vertices[i + 1].first*x_axis)/2.0;

			}
			else if (phase == 2) {
				vertical_area2.push_back(area());
				index = vertical_area2.size() - 1;
				vertical_area2[index].max_area = (2 * y_peak-vertices[i].second - vertices[i + 1].second)*(vertices[i].first - vertices[i+1].first)/2.0;
				vertical_area2[index].to = vertices[i+1];
				vertical_area2[index].from = vertices[i];
				vertical_area2[index].coeff_2 = gradient / 2.0;
				vertical_area2[index].coeff_1 = (constant + vertices[i].second - 2 * y_peak - gradient * vertices[i].first) / 2.0;
				vertical_area2[index].coeff_con = (2 * y_peak*vertices[i].first - constant * vertices[i].first - vertices[i].first*vertices[i].second)/2.0;

			}

		}
		phase = 1;
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].second - y_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i+1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					horizon_area1.push_back(area());
					index = horizon_area1.size() - 1;
					horizon_area1[index].from = vertices[i];
					horizon_area1[index].to = vertices[i + 1];
				}
				else {
					horizon_area2.push_back(area());
					index = horizon_area2.size() - 1;
					horizon_area2[index].to = vertices[i + 1];
					horizon_area2[index].from = vertices[i];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				horizon_area1.push_back(area());
				index = horizon_area1.size() - 1;
				horizon_area1[index].max_area = (2 * x_peak - vertices[i].first - vertices[i + 1].first)*(vertices[i + 1].second - vertices[i].second) / 2.0;
				horizon_area1[index].from = vertices[i];
				horizon_area1[index].to = vertices[i + 1];
				horizon_area1[index].coeff_2 = 1 / (2 * gradient);
				horizon_area1[index].coeff_1 = (vertices[i + 1].first - 2 * x_peak - (vertices[i + 1].second + constant) / gradient) / 2.0;
				horizon_area1[index].coeff_con = vertices[i + 1].second*(2 * x_peak + constant / gradient - vertices[i + 1].first)/2.0;

			}
			else if (phase == 2) {
				horizon_area2.push_back(area());
				index = horizon_area2.size() - 1;
				horizon_area2[index].max_area = (vertices[i].first - vertices[i + 1].first - 2 * y_axis)*(vertices[i].second - vertices[i + 1].second) / 2.0;
				horizon_area2[index].to = vertices[i + 1];
				horizon_area2[index].from = vertices[i];
				horizon_area2[index].coeff_2 = -1 / (2 * gradient);
				horizon_area2[index].coeff_1 = ((vertices[i].second + constant) / gradient - vertices[i].first+2*y_axis) / 2.0;
				horizon_area2[index].coeff_con = vertices[i].second*(vertices[i].first - constant / gradient - 2 * y_axis) / 2.0;
			}
		}
		
		
	}
	else if (quadrant == 2) {
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].first - x_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i+1].second-vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					vertical_area1.push_back(area());
					index = vertical_area1.size() - 1;
					vertical_area1[index].from = vertices[i];
					vertical_area1[index].to = vertices[i + 1];
				}
				else if (phase == 2) {
					vertical_area2.push_back(area());
					index = vertical_area2.size() - 1;
					vertical_area2[index].from = vertices[i];
					vertical_area2[index].to = vertices[i + 1];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				vertical_area1.push_back(area());
				index = vertical_area1.size() - 1;
				vertical_area1[index].from = vertices[i];
				vertical_area1[index].to = vertices[i + 1];
				vertical_area1[index].max_area = (2 * y_peak - vertices[i].second - vertices[i + 1].second)*(vertices[i].first - vertices[i + 1].first) / 2.0;
				vertical_area1[index].coeff_2 = -gradient / 2.0;
				vertical_area1[index].coeff_1 = (2 * y_peak -constant - vertices[i + 1].second + gradient * vertices[i + 1].first)/2.0;
				vertical_area1[index].coeff_con =  -vertices[i+1].first*(2*y_peak-constant-vertices[i+1].second)/ 2.0;

			}
			else if (phase == 2) {
				vertical_area2.push_back(area());
				index = vertical_area2.size() - 1;
				vertical_area2[index].from = vertices[i];
				vertical_area2[index].to = vertices[i + 1];
				vertical_area2[index].max_area = (vertices[i].second + vertices[i + 1].second - 2 * x_axis)*(vertices[i + 1].first - vertices[i].first)/2.0;
				vertical_area2[index].coeff_2 = gradient / 2.0;
				vertical_area2[index].coeff_1 = (constant + vertices[i].second - 2 * x_axis - gradient * vertices[i].first) / 2.0;
				vertical_area2[index].coeff_con = -vertices[i].first*(constant + vertices[i].second - 2 * x_axis)/2.0;
			}
		}
		phase = 1;
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].second - y_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i + 1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					horizon_area1.push_back(area());
					index = horizon_area1.size() - 1;
					horizon_area1[index].from = vertices[i];
					horizon_area1[index].to = vertices[i + 1];
				}
				else if (phase == 2) {
					horizon_area2.push_back(area());
					index = horizon_area2.size() - 1;
					horizon_area2[index].from = vertices[i];
					horizon_area2[index].to = vertices[i + 1];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				horizon_area1.push_back(area());
				index = horizon_area1.size() - 1;
				horizon_area1[index].from = vertices[i];
				horizon_area1[index].to = vertices[i + 1];
				horizon_area1[index].max_area = (2 * y_axis - vertices[i].first - vertices[i + 1].first)*(vertices[i + 1].second - vertices[i].second) / 2.0;
				horizon_area1[index].coeff_2 = 1 / (gradient*2.0);
				horizon_area1[index].coeff_1 = (-(constant + vertices[i + 1].second) / gradient - 2 * y_axis + vertices[i + 1].first) / 2.0;
				horizon_area1[index].coeff_con = vertices[i + 1].second*(2 * y_axis + constant / gradient - vertices[i + 1].first) / 2.0;

			}
			else if (phase == 2) {
				horizon_area2.push_back(area());
				index = horizon_area2.size() - 1;
				horizon_area2[index].from = vertices[i];
				horizon_area2[index].to = vertices[i + 1];
				horizon_area2[index].max_area = (vertices[i + 1].first + vertices[i].first - 2 * x_peak)*(vertices[i].second - vertices[i + 1].second) / 2.0;
				horizon_area2[index].coeff_2 = -1 / (gradient*2.0);
				horizon_area2[index].coeff_1 = (constant / gradient - vertices[i].first + 2 * x_peak + vertices[i].second / gradient) / 2.0;
				horizon_area2[index].coeff_con = vertices[i].second*(vertices[i].first - constant / gradient - 2 * x_peak) / 2.0;
			}
		}

			

	}
	else if (quadrant == 3) {
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].first - x_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i + 1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					vertical_area1.push_back(area());
					index = vertical_area1.size() - 1;
					vertical_area1[index].from = vertices[i];
					vertical_area1[index].to = vertices[i + 1];
				}
				else if (phase == 2) {
					vertical_area2.push_back(area());
					index = vertical_area2.size() - 1;
					vertical_area2[index].from = vertices[i];
					vertical_area2[index].to = vertices[i + 1];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				vertical_area1.push_back(area());
				index = vertical_area1.size() - 1;
				vertical_area1[index].from = vertices[i];
				vertical_area1[index].to = vertices[i + 1];
				vertical_area1[index].max_area = (2 * x_axis - vertices[i].second - vertices[i + 1].second)*(vertices[i].first - vertices[i + 1].first) / 2.0;
				vertical_area1[index].coeff_2 = -gradient / 2.0;
				vertical_area1[index].coeff_1 = (2*x_axis -constant -vertices[i+1].second+gradient*vertices[i+1].first) / 2.0;
				vertical_area1[index].coeff_con = -vertices[i + 1].first*(2*x_axis-constant-vertices[i+1].second) / 2.0;

			}
			else if (phase == 2) {
				vertical_area2.push_back(area());
				index = vertical_area2.size() - 1;
				vertical_area2[index].from = vertices[i];
				vertical_area2[index].to = vertices[i + 1];
				vertical_area2[index].max_area = (vertices[i].second + vertices[i + 1].second - 2 * y_peak)*(vertices[i + 1].first - vertices[i].first) / 2.0;
				vertical_area2[index].coeff_2 = gradient / 2.0;
				vertical_area2[index].coeff_1 = (constant + vertices[i].second - 2 * y_peak - gradient * vertices[i].first) / 2.0;
				vertical_area2[index].coeff_con = -vertices[i].first*(constant + vertices[i].second - 2 * y_peak) / 2.0;
			}
		}
		phase = 1;
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].second - y_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i + 1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					horizon_area1.push_back(area());
					index = horizon_area1.size() - 1;
					horizon_area1[index].from = vertices[i];
					horizon_area1[index].to = vertices[i + 1];
				}
				else if (phase == 2) {
					horizon_area2.push_back(area());
					index = horizon_area2.size() - 1;
					horizon_area2[index].from = vertices[i];
					horizon_area2[index].to = vertices[i + 1];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				horizon_area1.push_back(area());
				index = horizon_area1.size() - 1;
				horizon_area1[index].from = vertices[i];
				horizon_area1[index].to = vertices[i + 1];
				horizon_area1[index].max_area = (vertices[i + 1].first + vertices[i].first - 2 * x_peak)*(vertices[i].second - vertices[i+1].second)/2.0;
				horizon_area1[index].coeff_2 = 1 / (2 * gradient);
				horizon_area1[index].coeff_1 = (-(constant + vertices[i + 1].second) / gradient + vertices[i + 1].first - 2 * x_peak)/2.0;
				horizon_area1[index].coeff_con = -vertices[i + 1].second*(-constant / gradient + vertices[i + 1].first - 2 * x_peak)/2.0;
			}
			else if (phase == 2) {
				horizon_area2.push_back(area());
				index = horizon_area2.size() - 1;
				horizon_area2[index].from = vertices[i];
				horizon_area2[index].to = vertices[i + 1];
				horizon_area2[index].max_area = (y_axis * 2 - vertices[i].first - vertices[i + 1].first)*(vertices[i + 1].second - vertices[i].second) / 2.0;
				horizon_area2[index].coeff_2 = -1 / (2 * gradient);
				horizon_area2[index].coeff_1 = (2 * y_axis + (constant + vertices[i].second) / gradient - vertices[i].first) / 2.0;
				horizon_area2[index].coeff_con = -vertices[i].second*(2 * y_axis + constant / gradient - vertices[i].first)/2.0;
			}
		}

	}
	else {
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].first - x_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i + 1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					vertical_area1.push_back(area());
					index = vertical_area1.size() - 1;
					vertical_area1[index].from = vertices[i];
					vertical_area1[index].to = vertices[i + 1];
				}
				else {
					vertical_area2.push_back(area());
					index = vertical_area2.size() - 1;
					vertical_area2[index].to = vertices[i + 1];
					vertical_area2[index].from = vertices[i];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				vertical_area1.push_back(area());
				index = vertical_area1.size() - 1;
				vertical_area1[index].max_area = (vertices[i + 1].second + vertices[i].second - 2 * y_peak)*(vertices[i + 1].first - vertices[i].first)/2.0;
				vertical_area1[index].from = vertices[i];
				vertical_area1[index].to = vertices[i + 1];
				vertical_area1[index].coeff_2 = -gradient / 2.0;
				vertical_area1[index].coeff_1 = (gradient*vertices[i + 1].first - constant - vertices[i + 1].second + 2 * y_peak)/2.0;
				vertical_area1[index].coeff_con = vertices[i + 1].first*(constant + vertices[i + 1].second - 2 * y_peak)/2.0;
			}
			else if (phase == 2) {
				vertical_area2.push_back(area());
				index = vertical_area2.size() - 1;
				vertical_area2[index].from = vertices[i];
				vertical_area2[index].to = vertices[i + 1];
				vertical_area2[index].max_area = (2 * x_axis - vertices[i + 1].second - vertices[i].second)*(vertices[i].first - vertices[i + 1].first) / 2.0;
				vertical_area2[index].coeff_2 = gradient / 2.0;
				vertical_area2[index].coeff_1 = (-gradient * vertices[i].first - 2 * x_axis + constant + vertices[i].second)/2.0;
				vertical_area2[index].coeff_con = vertices[i].first*(2 * x_axis - constant - vertices[i].second) / 2.0;
			}
		}
		phase = 1;
		for (int i = 0; i < vertices.size() - 1; i++) {
			if (abs(vertices[i].second - y_peak)<=DBL_EPSILON)
				phase = 2;
			if (abs(vertices[i + 1].first - vertices[i].first)<=DBL_EPSILON || abs(vertices[i + 1].second - vertices[i].second)<=DBL_EPSILON) {
				if (phase == 1) {
					horizon_area1.push_back(area());
					index = horizon_area1.size() - 1;
					horizon_area1[index].from = vertices[i];
					horizon_area1[index].to = vertices[i + 1];
				}
				else {
					horizon_area2.push_back(area());
					index = horizon_area2.size() - 1;
					horizon_area2[index].to = vertices[i + 1];
					horizon_area2[index].from = vertices[i];
				}
				continue;
			}
			gradient = (vertices[i + 1].second - vertices[i].second) / (vertices[i + 1].first - vertices[i].first);
			constant = vertices[i].second - gradient * vertices[i].first;
			if (phase == 1) {
				horizon_area1.push_back(area());
				index = horizon_area1.size() - 1;
				horizon_area1[index].from = vertices[i];
				horizon_area1[index].to = vertices[i + 1];
				horizon_area1[index].max_area = (vertices[i + 1].first + vertices[i].first - 2 * y_axis)*(vertices[i].second - vertices[i + 1].second) / 2.0;
				horizon_area1[index].coeff_2 = 1 / (gradient * 2);
				horizon_area1[index].coeff_1 = (-(constant + vertices[i + 1].second) / gradient + vertices[i + 1].first - 2 * y_axis)/2.0;
				horizon_area1[index].coeff_con = -vertices[i + 1].second*(-constant / gradient + vertices[i + 1].first - 2 * y_axis)/2.0;
			}
			else if (phase == 2) {
				horizon_area2.push_back(area());
				index = horizon_area2.size() - 1;
				horizon_area2[index].from = vertices[i];
				horizon_area2[index].to = vertices[i + 1];
				horizon_area2[index].max_area = (2 * x_peak - vertices[i + 1].first - vertices[i].first)*(vertices[i + 1].second - vertices[i].second) / 2.0;
				horizon_area2[index].coeff_2 = -1 / (2 * gradient);
				horizon_area2[index].coeff_1 = ((constant + vertices[i].second) / gradient + 2 * x_peak - vertices[i].first) / 2.0;
				horizon_area2[index].coeff_con = -vertices[i].second*(2 * x_peak + constant / gradient - vertices[i].first)/2.0;
			}
		}

	}

	for (int i = vertical_area1.size() - 2; i >= 0; i--) {
		vertical_area1[i].accumulated_area = vertical_area1[i + 1].max_area + vertical_area1[i + 1].accumulated_area;
	}
	for (int i = 1; i < vertical_area2.size(); i++) {
		vertical_area2[i].accumulated_area = vertical_area2[i - 1].max_area + vertical_area2[i - 1].accumulated_area;
	}
	for (int i = horizon_area1.size() - 2; i >= 0; i--) {
		horizon_area1[i].accumulated_area = horizon_area1[i + 1].max_area + horizon_area1[i + 1].accumulated_area;
	}
	for (int i = 1; i < horizon_area2.size(); i++) {
		horizon_area2[i].accumulated_area = horizon_area2[i - 1].max_area + horizon_area2[i - 1].accumulated_area;
	}
	
	empty = 0;

	for (int i = 0; i < vertical_area1.size(); i++)
		empty += vertical_area1[i].max_area;
	for (int i = 0; i < vertical_area2.size(); i++)
		empty += vertical_area2[i].max_area;
	
}

double area_struct::cal_area(double input_x, double input_y, int v1_idx, int v2_idx, int h1_idx , int h2_idx ) {
	
	double area = 0;
	if (input_x < y_axis && input_x < x_peak || input_x>y_axis && input_x>x_peak)
		input_x = x_peak;
	if (input_y<x_axis && input_y < y_peak || input_y>x_axis && input_y>y_peak)
		input_y = y_peak;
	
	if (v1_idx == -2) {
		v1_idx = find_index(input_x, 1);
	}
	
	if (v2_idx == -2) {
		v2_idx = find_index(input_x, 2);
	}

	if (h1_idx == -2) {
		h1_idx = find_index(input_y, 3);
	}

	if (h2_idx == -2) {
		h2_idx = find_index(input_y, 4);
	}

	if (quadrant == 1) {

		if (check_external(vertical_area2[v2_idx].from, vertical_area2[v2_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				area -= abs((y_peak - input_y)*(x_peak - input_x));
				return abs((input_y - x_axis)*(input_x - y_axis)) - (empty - area);
				/*
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				area -= (y_peak - input_y)*(x_peak - input_x);
				return (input_y - x_axis)*(input_x - y_axis) - (empty - area);
				*/
			}
			else {
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}
				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				return abs((input_x - y_axis)*(y_peak - x_axis)) - (empty - area);
				/*
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				return (input_x - y_axis)*(y_peak - x_axis) - (empty - area);*/
			}

		}
		else {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {
				
				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				
				return abs((input_y - x_axis)*(x_peak - y_axis)) - (empty - area);
				/*
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				return (input_y - x_axis)*(x_peak - y_axis) - (empty - area);*/
			}
			else {
				area = abs((input_x - y_axis)*(input_y - x_axis));
				if (v1_idx != -1) {
					area -= (vertical_area1[0].accumulated_area + vertical_area1[0].max_area) - (vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x));
				}
				if (h2_idx != -1) {
					area -= (horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area) - (horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y));
				}
				return area;
			}
		}
		

	}
	else if (quadrant == 2) {
		if (check_external(vertical_area1[v1_idx].from, vertical_area1[v1_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				area -= abs((y_peak - input_y)*(x_peak - input_x));
				return abs((input_y - x_axis)*(input_x - y_axis)) - (empty - area);

			}
			else {
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}
				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if(vertical_area2.size()!=0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}
				
				return abs((input_x - y_axis)*(y_peak - x_axis)) - (empty - area);
			}
		}
		else {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1){
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				
				return abs((input_y - x_axis)*(x_peak - y_axis)) - (empty - area);
			}
			else {
				area = abs((y_axis - input_x)*(input_y - x_axis));
				if (v2_idx != -1) {
					area -= (vertical_area2[vertical_area2.size()-1].accumulated_area + vertical_area2[vertical_area2.size()-1].max_area) - (vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x));
				}
				if (h1_idx != -1) {
					area -= (horizon_area1[0].accumulated_area + horizon_area1[0].max_area) - (horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y));
				}
				return area;
			}
		}
	}
	else if (quadrant == 3) {
		if (check_external(vertical_area2[v2_idx].from, vertical_area2[v2_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				area -= abs((y_peak - input_y)*(x_peak - input_x));
				return abs((input_y - x_axis)*(input_x - y_axis)) - (empty - area);

			}
			else {
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}
				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				return abs((input_x - y_axis)*(y_peak - x_axis)) - (empty - area);
			}
		}
		else {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				return abs((input_y - x_axis)*(x_peak - y_axis)) - (empty - area);
			}
			else {
				area = abs((y_axis - input_x)*(input_y - x_axis));
				if (v1_idx != -1) {
					area -= (vertical_area1[0].accumulated_area + vertical_area1[0].max_area) - (vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x));
				}
				if (h2_idx != -1) {
					area -= (horizon_area2[horizon_area2.size()-1].accumulated_area + horizon_area2[horizon_area2.size()-1].max_area) - (horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y));
				}
				return area;
			}
		}
	}
	else if (quadrant == 4) {
		if (check_external(vertical_area1[v1_idx].from, vertical_area1[v1_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				area -= abs((y_peak - input_y)*(x_peak - input_x));
				return abs((input_y - x_axis)*(input_x - y_axis)) - (empty - area);

			}
			else {
				if (v1_idx != -1) {
					area += vertical_area1[v1_idx].accumulated_area + vertical_area1[v1_idx].partial_area(input_x);
				}
				else {
					if (vertical_area1.size() != 0)
						area += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}
				if (v2_idx != -1) {
					area += vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x);
				}
				else {
					if (vertical_area2.size() != 0)
						area += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				return abs((input_x - y_axis)*(y_peak - x_axis)) - (empty - area);
			}
		}
		else {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1) {
					area += horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y);
				}
				else {
					if (horizon_area1.size() != 0) {
						area += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}
				if (h2_idx != -1) {
					area += horizon_area2[h2_idx].accumulated_area + horizon_area2[h2_idx].partial_area(input_y);
				}
				else {
					if (horizon_area2.size() != 0) {
						area += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				return abs((input_y - x_axis)*(x_peak - y_axis)) - (empty - area);
			}
			else {
				area = abs((y_axis - input_x)*(input_y - x_axis));
				if (v2_idx != -1) {
					area -= (vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area) - (vertical_area2[v2_idx].accumulated_area + vertical_area2[v2_idx].partial_area(input_x));
				}
				if (h1_idx != -1) {
					area -= (horizon_area1[0].accumulated_area + horizon_area1[0].max_area) - (horizon_area1[h1_idx].accumulated_area + horizon_area1[h1_idx].partial_area(input_y));
				}
				return area;
			}
		}
	}
	
}

int area_struct::find_index(double value, int type) { // type 1: vertical_area1 / type 2: vertical_area2 / type 3: horizon_area1 / type 4: horizon_area2
	
	if (type < 3) {
		if (value< y_axis && value < x_peak || value>y_axis &&value>x_peak)
			value = x_peak;
	}
	else {
		if (value<x_axis && value < y_peak || value>x_axis && value>y_peak)
			value = y_peak;
	}
	
	
	if (type == 1) {
		if (vertical_area1.size() == 0)
			return -1;
		if (quadrant == 1 || quadrant == 4)
			return binary_search_index(vertical_area1, value, 1, 0, vertical_area1.size() - 1, 1);
		else
			return binary_search_index(vertical_area1, value, 1, 0, vertical_area1.size() - 1, 2);
	}
	else if (type == 2) {
		if (vertical_area2.size() == 0)
			return -1;
		if (quadrant == 2 || quadrant == 3)
			return binary_search_index(vertical_area2, value, 1, 0, vertical_area2.size() - 1, 1);
		else
			return binary_search_index(vertical_area2, value, 1, 0, vertical_area2.size() - 1, 2);
	}
	else if (type == 3) {
		if (horizon_area1.size() == 0)
			return -1;
		if (quadrant == 1 || quadrant == 2)
			return binary_search_index(horizon_area1, value, 2, 0, horizon_area1.size() - 1, 1);
		else
			return binary_search_index(horizon_area1, value, 2, 0, horizon_area1.size() - 1, 2);
	}
	else if (type == 4) {
		if (horizon_area2.size() == 0)
			return -1;
		if (quadrant == 3 || quadrant == 4)
			return binary_search_index(horizon_area2, value, 2, 0, horizon_area2.size() - 1, 1);
		else
			return binary_search_index(horizon_area2, value, 2, 0, horizon_area2.size() - 1, 2);
	}
}

int binary_search_index(vector<area>& vec, double value, int type,int start, int end, int order) {// type 1: vertical / type 2: horizon // order 1: increasing / order 2: decreasing
	if (end < start) //can not find
		return -1;

	int mid = (start + end) / 2;

	if (type == 1) {
		if (vec[mid].from.first <= value && vec[mid].to.first >= value || vec[mid].from.first >= value && vec[mid].to.first <= value)
			return mid;
		else if (max(vec[mid].from.first, vec[mid].to.first) < value) {
			if (order == 1)
				return binary_search_index(vec, value, type, mid + 1, end, order);
			else if (order == 2)
				return binary_search_index(vec, value, type, start, mid - 1, order);
		}
		else {
			if (order == 1)
				return binary_search_index(vec, value, type, start, mid - 1, order);
			else if (order == 2)
				return binary_search_index(vec, value, type, mid + 1, end, order);
		}
	}
	else if (type == 2) {
		if (max(vec[mid].from.second, vec[mid].to.second) >= value && min(vec[mid].from.second, vec[mid].to.second) <= value)
			return mid;
		else if (max(vec[mid].from.second, vec[mid].to.second) < value)
		{
			if (order == 1)
				return binary_search_index(vec, value, type, mid + 1, end, order);
			else if (order == 2)
				return binary_search_index(vec, value, type, start, mid - 1, order);
		}
		else{
			if (order == 1)
				return binary_search_index(vec, value, type, start, mid - 1, order);
			else if (order == 2)
				return binary_search_index(vec, value, type, mid + 1, end, order);

		}
	}
	
}

void area_struct::cal_coeff(double input_x, double input_y, double* coeff_x2, double* coeff_x1, double* coeff_y2, double* coeff_y1, double* coeff_xy,double* coeff_con) {
	int v1_idx, v2_idx, h1_idx, h2_idx;
	bool x_check = false;
	bool y_check = false;
	//double area = 0;
	if (input_x < y_axis && input_x < x_peak || input_x>y_axis && input_x>x_peak)
	{
		input_x = x_peak;
		x_check = true;
	}
	if (input_y<x_axis && input_y < y_peak || input_y>x_axis && input_y>y_peak)
	{
		input_y = y_peak;
		y_check = true;
	}

	v1_idx = find_index(input_x, 1);
	v2_idx = find_index(input_x, 2);
	h1_idx = find_index(input_y, 3);
	h2_idx = find_index(input_y, 4);
	
	if (quadrant == 1) {

		if (check_external(vertical_area2[v2_idx].from, vertical_area2[v2_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {
				
				if (v1_idx != -1) {//general version
					if(!x_check) {
						(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
						(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
						(*coeff_con) += vertical_area1[v1_idx].coeff_con;
						(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
					}

					
				}
				else {
					if(vertical_area1.size()!=0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					if(!x_check) {
						(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
						(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
						(*coeff_con) += vertical_area2[v2_idx].coeff_con;
						(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
					}
			
				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					if(!y_check) {
						(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
						(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
						(*coeff_con) += horizon_area1[h1_idx].coeff_con;
						(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
					}
					
				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					if(!y_check) {
						(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
						(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
						(*coeff_con) += horizon_area2[h2_idx].coeff_con;
						(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
					}
					
				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				if (x_check) {
					(*coeff_con) += x_peak * (y_peak - x_axis);
				}
				else {
					(*coeff_x1) += (y_peak - x_axis);
				}

				if (y_check) {
					(*coeff_con) += y_peak * (x_peak - y_axis);
				}
				else {
					(*coeff_y1) += (x_peak - y_axis);
				}
				
				(*coeff_con) += (y_axis*x_axis - x_peak * y_peak - empty);
				

			}
			else {
				if (v1_idx != -1) {
					
						(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
						(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
						(*coeff_con) += vertical_area1[v1_idx].coeff_con;
						(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
					
				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					
						(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
						(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
						(*coeff_con) += vertical_area2[v2_idx].coeff_con;
						(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
					
				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}
				
				(*coeff_x1) += (y_peak - x_axis);
				(*coeff_con) += -y_axis*(y_peak - x_axis) - empty;
			}

		}
		else {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {

				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				(*coeff_y1) += (x_peak - y_axis);
				(*coeff_con) += -x_axis * (x_peak - y_axis) - empty;
			}
			else {

				if (v1_idx != -1) {
					(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
					(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
					(*coeff_con) += vertical_area1[v1_idx].coeff_con;
					(*coeff_con) -= (vertical_area1[0].accumulated_area + vertical_area1[0].max_area) - (vertical_area1[v1_idx].accumulated_area);

				}
				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) -= (horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area) - (horizon_area2[h2_idx].accumulated_area );

				}
				(*coeff_xy) += 1;
				(*coeff_x1) -= x_axis;
				(*coeff_y1) -= y_axis;
				(*coeff_con) += x_axis * y_axis;
			}
		}


	}
	else if (quadrant == 2) {
		if (check_external(vertical_area1[v1_idx].from, vertical_area1[v1_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					if (!x_check) {
						(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
						(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
						(*coeff_con) += vertical_area1[v1_idx].coeff_con;
						(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
					}


				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					if (!x_check) {
						(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
						(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
						(*coeff_con) += vertical_area2[v2_idx].coeff_con;
						(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
					}

				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
						(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
						(*coeff_con) += horizon_area1[h1_idx].coeff_con;
						(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
						(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
						(*coeff_con) += horizon_area2[h2_idx].coeff_con;
						(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				if (x_check) {
					(*coeff_con) += x_peak * (-y_peak + x_axis);
				}
				else {
					(*coeff_x1) += (-y_peak + x_axis);
				}
				if (y_check) {
					(*coeff_con) += y_peak*(y_axis - x_peak);
				}
				else {
					(*coeff_y1) += (-x_peak + y_axis);
				}
				(*coeff_con) += (x_peak*y_peak - y_axis * x_axis - empty);

			}
			else {
				if (v1_idx != -1) {
					(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
					(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
					(*coeff_con) += vertical_area1[v1_idx].coeff_con;
					(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
					(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
					(*coeff_con) += vertical_area2[v2_idx].coeff_con;
					(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}
				(*coeff_x1) -= (y_peak - x_axis);
				(*coeff_con) += y_axis * (y_peak - x_axis) - empty;

			}
		}
		else {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				(*coeff_y1) -= (x_peak - y_axis);
				(*coeff_con) += -x_axis * (y_axis - x_peak) - empty;

			}
			else {
				if (v2_idx != -1) {
					(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
					(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
					(*coeff_con) += vertical_area2[v2_idx].coeff_con;
					(*coeff_con) -= (vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area) - (vertical_area2[v2_idx].accumulated_area);

				}
				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) -= (horizon_area1[0].accumulated_area + horizon_area1[0].max_area) - (horizon_area1[h1_idx].accumulated_area);

				}
				(*coeff_xy) -= 1;
				(*coeff_x1) += x_axis;
				(*coeff_y1) += y_axis;
				(*coeff_con) -= x_axis * y_axis;
			}
		}
	}
	else if (quadrant == 3) {
		if (check_external(vertical_area2[v2_idx].from, vertical_area2[v2_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {
				if (v1_idx != -1) {//general version
					if (!x_check) {
						(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
						(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
						(*coeff_con) += vertical_area1[v1_idx].coeff_con;
						(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
					}


				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					if (!x_check) {
						(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
						(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
						(*coeff_con) += vertical_area2[v2_idx].coeff_con;
						(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
					}

				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
						(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
						(*coeff_con) += horizon_area1[h1_idx].coeff_con;
						(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
						(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
						(*coeff_con) += horizon_area2[h2_idx].coeff_con;
						(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				if (x_check) {
					(*coeff_con) += x_peak * (y_peak - x_axis);
				}
				else {
					(*coeff_x1) += (y_peak - x_axis);
				}
				if (y_check) {
					(*coeff_con) += y_peak * (x_peak - y_axis);
				
				}
				else {
					(*coeff_y1) += (x_peak - y_axis);
				}
				(*coeff_con) += (y_axis*x_axis - x_peak * y_peak - empty);

			}
			else {
				if (v1_idx != -1) {
					(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
					(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
					(*coeff_con) += vertical_area1[v1_idx].coeff_con;
					(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
					(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
					(*coeff_con) += vertical_area2[v2_idx].coeff_con;
					(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}
				(*coeff_x1) += (y_peak - x_axis);
				(*coeff_con) += y_axis * (x_axis - y_peak) - empty;
			}
		}
		else {
			if (check_external(horizon_area1[h1_idx].from, horizon_area1[h1_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				(*coeff_y1) += (x_peak - y_axis);
				(*coeff_con) += x_axis * (y_axis - x_peak) - empty;

			}
			else {
				
				if (v1_idx != -1) {
					(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
					(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
					(*coeff_con) += vertical_area1[v1_idx].coeff_con;
					(*coeff_con) -= (vertical_area1[0].accumulated_area + vertical_area1[0].max_area) - (vertical_area1[v1_idx].accumulated_area );

				}
				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) -= (horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area) - (horizon_area2[h2_idx].accumulated_area);


				}
				(*coeff_xy) += 1;
				(*coeff_x1) -= x_axis;
				(*coeff_y1) -= y_axis;
				(*coeff_con) += x_axis * y_axis;
			}
		}
	}
	else if (quadrant == 4) {
		if (check_external(vertical_area1[v1_idx].from, vertical_area1[v1_idx].to, pair<double, double>(input_x, input_y))) {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {

				if (v1_idx != -1) {//general version
					if (!x_check) {
						(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
						(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
						(*coeff_con) += vertical_area1[v1_idx].coeff_con;
						(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
					}


				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					if (!x_check) {
						(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
						(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
						(*coeff_con) += vertical_area2[v2_idx].coeff_con;
						(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
					}

				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}

				if (h1_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
						(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
						(*coeff_con) += horizon_area1[h1_idx].coeff_con;
						(*coeff_con) += horizon_area1[h1_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					if (!y_check) {
						(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
						(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
						(*coeff_con) += horizon_area2[h2_idx].coeff_con;
						(*coeff_con) += horizon_area2[h2_idx].accumulated_area;
					}

				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}

				if (x_check) {
					(*coeff_con) += x_peak * (x_axis - y_peak);
				}
				else
				{
					(*coeff_x1) += (-y_peak + x_axis);
				}
				if (y_check) {
					(*coeff_con) += y_peak * (y_axis - x_peak);
				}
				else {
					(*coeff_y1) += (-x_peak + y_axis);
				}
				
				(*coeff_con) += (-x_axis * y_axis + x_peak * y_peak - empty);

			}
			else {
				if (v1_idx != -1) {
					(*coeff_x2) += vertical_area1[v1_idx].coeff_2;
					(*coeff_x1) += vertical_area1[v1_idx].coeff_1;
					(*coeff_con) += vertical_area1[v1_idx].coeff_con;
					(*coeff_con) += vertical_area1[v1_idx].accumulated_area;
				}
				else {
					if (vertical_area1.size() != 0)
						(*coeff_con) += max(vertical_area1[0].accumulated_area + vertical_area1[0].max_area, vertical_area1[vertical_area1.size() - 1].accumulated_area + vertical_area1[vertical_area1.size() - 1].max_area);
				}

				if (v2_idx != -1) {
					(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
					(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
					(*coeff_con) += vertical_area2[v2_idx].coeff_con;
					(*coeff_con) += vertical_area2[v2_idx].accumulated_area;
				}
				else {
					if (vertical_area2.size() != 0)
						(*coeff_con) += max(vertical_area2[0].accumulated_area + vertical_area2[0].max_area, vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area);
				}
				(*coeff_x1) -= (y_peak - x_axis);
				(*coeff_con) += -y_axis * (x_axis - y_peak) - empty;
			}
		}
		else {
			if (check_external(horizon_area2[h2_idx].from, horizon_area2[h2_idx].to, pair<double, double>(input_x, input_y))) {
				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
				}
				else {
					if (horizon_area1.size() != 0) {
						(*coeff_con) += max(horizon_area1[0].accumulated_area + horizon_area1[0].max_area, horizon_area1[horizon_area1.size() - 1].accumulated_area + horizon_area1[horizon_area1.size() - 1].max_area);
					}
				}

				if (h2_idx != -1) {
					(*coeff_y2) += horizon_area2[h2_idx].coeff_2;
					(*coeff_y1) += horizon_area2[h2_idx].coeff_1;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
					(*coeff_con) += horizon_area2[h2_idx].coeff_con;
				}
				else {
					if (horizon_area2.size() != 0) {
						(*coeff_con) += max(horizon_area2[0].accumulated_area + horizon_area2[0].max_area, horizon_area2[horizon_area2.size() - 1].accumulated_area + horizon_area2[horizon_area2.size() - 1].max_area);
					}
				}
				(*coeff_y1) -= (x_peak - y_axis);
				(*coeff_con) += x_axis * (x_peak - y_axis) - empty;
			}
			else {
				if (v2_idx != -1) {
					(*coeff_x2) += vertical_area2[v2_idx].coeff_2;
					(*coeff_x1) += vertical_area2[v2_idx].coeff_1;
					(*coeff_con) += vertical_area2[v2_idx].coeff_con;
					(*coeff_con) -= (vertical_area2[vertical_area2.size() - 1].accumulated_area + vertical_area2[vertical_area2.size() - 1].max_area) - (vertical_area2[v2_idx].accumulated_area);

				}
				if (h1_idx != -1) {
					(*coeff_y2) += horizon_area1[h1_idx].coeff_2;
					(*coeff_y1) += horizon_area1[h1_idx].coeff_1;
					(*coeff_con) += horizon_area1[h1_idx].coeff_con;
					(*coeff_con) -= (horizon_area1[0].accumulated_area + horizon_area1[0].max_area) - (horizon_area1[h1_idx].accumulated_area);

				}
				(*coeff_xy) -= 1;
				(*coeff_x1) += x_axis;
				(*coeff_y1) += y_axis;
				(*coeff_con) -= x_axis * y_axis;
			}
		}
	}

}

bool check_external(pair<double, double> from, pair<double, double> to, pair<double, double> point) {
	pair<double, double> v1 = pair<double, double>(point.first - from.first, point.second - from.second);
	pair<double, double> v2 = pair<double, double>(to.first - from.first, to.second - from.second);
	return (v1.first*v2.second - v1.second*v2.first) > DBL_EPSILON;
}

void find_max_point(double x_low, double x_high, double y_low, double y_high, vertex* top, vertex* bot, vertex* left, vertex* right) {
	edge* e_list[8] = { NULL };//the list including the edges that possibly intersect with rectangle. from the left bottom vertex, ccw
	vector<int> e_check_list;// 0: not yet decided, 1: include, 2: not include 3: not intersect
	for (int i = 0; i < 8; i++) {
		e_check_list.push_back(0);
	}
	vertex* temp = bot;
	vector<pair<double, double>> domain;
	edge* e_temp = bot->edge_ccw;
	int phase = 1;
	int d_num = 0;
	pair<bool, pair<double, double>> p_temp;

	do {
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
					if (phase >= 3)
						e_list[0] = e_temp;
					else
						e_list[1] = e_temp;
				}

			}
			p_temp = e_temp->y_check((y_low + y_high) / 2 + height);
			if (p_temp.second.first >= x_high + width || p_temp.second.first <= x_low)
				p_temp.first = false;
			if (p_temp.first) {
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
					if (abs(e_temp->from->v_y - top->v_y) <= FLT_EPSILON) {
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
	double x_mid=x_low, y_mid=y_low;
	if (quad1.vertical_area2[quad1.vertical_area2.size() - 1].to.second < quad1.y_peak) {
		y_mid = quad1.vertical_area2[quad1.vertical_area2.size() - 1].to.second - height;
	}
	else {
		y_mid = quad4.vertical_area1[0].from.second;
	}

	if (quad1.horizon_area1[0].from.first < quad1.x_peak)
		x_mid = quad1.horizon_area1[0].from.first - width;
	else
		x_mid = quad3.horizon_area1[0].from.first;
	if (x_low<x_mid && x_mid<x_high && abs(x_low - x_mid)>DBL_EPSILON && abs(x_mid - x_high)>DBL_EPSILON) {
		if (y_low<y_mid && y_mid<y_high &&abs(y_low - y_mid)>DBL_EPSILON && abs(y_mid - y_high)>DBL_EPSILON) {
			domain.push_back(pair<double, double>(x_low, y_low));
			domain.push_back(pair<double, double>(x_mid, y_low));
			domain.push_back(pair<double, double>(x_mid, y_mid));
			domain.push_back(pair<double, double>(x_low, y_mid));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
			domain.clear();

			domain.push_back(pair<double, double>(x_mid, y_low));
			domain.push_back(pair<double, double>(x_high, y_low));
			domain.push_back(pair<double, double>(x_high, y_mid));
			domain.push_back(pair<double, double>(x_mid, y_mid));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
			domain.clear();

			domain.push_back(pair<double, double>(x_low, y_mid));
			domain.push_back(pair<double, double>(x_mid, y_mid));
			domain.push_back(pair<double, double>(x_mid, y_high));
			domain.push_back(pair<double, double>(x_low, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
			domain.clear();

			domain.push_back(pair<double, double>(x_mid, y_mid));
			domain.push_back(pair<double, double>(x_high, y_mid));
			domain.push_back(pair<double, double>(x_high, y_high));
			domain.push_back(pair<double, double>(x_mid, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
		}
		else {
			domain.push_back(pair<double, double>(x_low, y_low));
			domain.push_back(pair<double, double>(x_mid, y_low));
			domain.push_back(pair<double, double>(x_mid, y_high));
			domain.push_back(pair<double, double>(x_low, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
			domain.clear();

			domain.push_back(pair<double, double>(x_mid, y_low));
			domain.push_back(pair<double, double>(x_high, y_low));
			domain.push_back(pair<double, double>(x_high, y_high));
			domain.push_back(pair<double, double>(x_mid, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
		}
	}
	else {
		if (y_low<y_mid && y_mid<y_high &&abs(y_low - y_mid)>DBL_EPSILON && abs(y_mid - y_high)>DBL_EPSILON) {
			domain.push_back(pair<double, double>(x_low, y_low));
			domain.push_back(pair<double, double>(x_high, y_low));
			domain.push_back(pair<double, double>(x_high, y_mid));
			domain.push_back(pair<double, double>(x_low, y_mid));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
			domain.clear();

			domain.push_back(pair<double, double>(x_low, y_mid));
			domain.push_back(pair<double, double>(x_high, y_mid));
			domain.push_back(pair<double, double>(x_high, y_high));
			domain.push_back(pair<double, double>(x_low, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
		}
		else {
			domain.push_back(pair<double, double>(x_low, y_low));
			domain.push_back(pair<double, double>(x_high, y_low));
			domain.push_back(pair<double, double>(x_high, y_high));
			domain.push_back(pair<double, double>(x_low, y_high));
			recursive_find(x_low, x_high, y_low, y_high, width, height, e_check_list, e_list, d_num, domain);
		}
	}
	
}

int edge_check(edge* e_temp, double left, double right, double top, double bot) {
	if ((e_temp->from->v_x <= left && e_temp->to->v_x <= left) || (e_temp->from->v_x >= right && e_temp->to->v_x >= right) || (e_temp->from->v_y <= bot && e_temp->to->v_y <= bot) || (e_temp->from->v_y >= top && e_temp->to->v_y >= top))
		return 0;
	else
		return 1;
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
		
		//edge* edge_temp = e_list[inter_vector[0]];
		coeff_set co_set1,co_set2,co_set3,co_set4;
		co_set1.quad = 1;
		co_set2.quad = 2;
		co_set3.quad = 3;
		co_set4.quad = 4;
		co_set1.x_peak = quad1.x_peak;
		co_set1.y_peak = quad1.y_peak;
		co_set2.x_peak = quad2.x_peak;
		co_set2.y_peak = quad2.y_peak;
		co_set3.x_peak = quad3.x_peak;
		co_set3.y_peak = quad3.y_peak;
		co_set4.x_peak = quad4.x_peak;
		co_set4.y_peak = quad4.y_peak;
		
		double c_x = 0.0;
		double c_y = 0.0;
		double c_xy = 0.0;
		double c_x_2 = 0.0;
		double c_y_2 = 0.0;
		double constant = 0.0;
		
		double sol_x;
		double sol_y;
		double temp_max_value;
		int t;
		pair<double, double> temp_pos = make_point(domain);
		
		quad1.cal_coeff(temp_pos.first+width, temp_pos.second+height, &co_set1.x2, &co_set1.x1, &co_set1.y2, &co_set1.y1, &co_set1.xy, &co_set1.constant);
		quad2.cal_coeff(temp_pos.first, temp_pos.second+height, &co_set2.x2, &co_set2.x1, &co_set2.y2, &co_set2.y1, &co_set2.xy, &co_set2.constant);
		quad3.cal_coeff(temp_pos.first, temp_pos.second, &co_set3.x2, &co_set3.x1, &co_set3.y2, &co_set3.y1, &co_set3.xy, &co_set3.constant);
		quad4.cal_coeff(temp_pos.first+width, temp_pos.second, &co_set4.x2, &co_set4.x1, &co_set4.y2, &co_set4.y1, &co_set4.xy, &co_set4.constant);

		c_x_2 += co_set1.x2 + co_set4.x2 + co_set2.x2 +co_set3.x2;
		c_y_2 += co_set2.y2 + co_set3.y2 + co_set1.y2 +co_set4.y2;
		c_x += 2 * width*co_set1.x2 + 2* width* co_set4.x2 + co_set1.x1 +co_set4.x1 + co_set2.x1+co_set3.x1+co_set1.xy*height+co_set2.xy*height;
		c_y += 2 * height*co_set2.y2 + 2 * height*co_set1.y2 + co_set2.y1 + co_set3.y1+co_set1.y1+co_set4.y1+co_set1.xy*width+co_set4.xy*width;
		c_xy += co_set1.xy + co_set2.xy + co_set3.xy + co_set4.xy;
		constant += co_set1.x2*width*width + co_set4.x2*width*width + co_set1.x1*width + co_set4.x1*width + co_set2.y2*height*height + co_set1.y2*height*height + co_set2.y1*height+co_set1.y1*height;
		constant += co_set1.constant + co_set2.constant + co_set3.constant + co_set4.constant + co_set1.xy*width*height;
		
		/*
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
			if (4 * c_x_2*c_y_2 - c_xy * c_xy > 0 && c_x_2<0) {
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
			i_temp = line_intersect(domain[k - 1].first, domain[k - 1].second, domain[k].first, domain[k].second, e_temp.temp_x1, e_temp.temp_y1, e_temp.temp_x2, e_temp.temp_y2);

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
		i_temp = line_intersect(domain[domain.size() - 1].first, domain[domain.size() - 1].second, domain[0].first, domain[0].second, e_temp.temp_x1, e_temp.temp_y1, e_temp.temp_x2, e_temp.temp_y2);
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
			i_temp = line_intersect(v_list2[k - 1].first, v_list2[k - 1].second, v_list2[k].first, v_list2[k].second, e_temp2.temp_x1, e_temp2.temp_y1, e_temp2.temp_x2, e_temp2.temp_y2);
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
		i_temp = line_intersect(v_list2[v_list2.size() - 1].first, v_list2[v_list2.size() - 1].second, v_list2[0].first, v_list2[0].second, e_temp2.temp_x1, e_temp2.temp_y1, e_temp2.temp_x2, e_temp2.temp_y2);
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

pair<double, double> line_intersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
	double a, b, c, d;
	if (abs(x1 - x2) <= DBL_EPSILON) {
		c = (y4 - y3) / (x4 - x3);
		return pair<double, double>(x1, c*(x1 - x3) + y3);
	}
	else if (abs(x3 - x4) <= DBL_EPSILON) {
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

double cal_area_coeff(double input_x, double input_y, double c_x2, double c_xy, double c_y2, double c_x1, double c_y1, double c_con) {
	return c_x2 * input_x*input_x + c_xy * input_x*input_y + c_y2 * input_y*input_y + c_x1 * input_x + c_y1 * input_y + c_con;
}

pair<double, double> make_point(vector<pair<double, double>> domain) {
	pair<double, double> temp;
	temp.first = (domain[0].first + domain[1].first) / 2.0;
	temp.second = (domain[0].second + domain[1].second) / 2.0;
	temp.first = (domain[2].first + temp.first) / 2.0;
	temp.second = (domain[2].second + temp.second) / 2.0;
	return temp;
}

void find_max_area_side(double x_2, double xy, double y_2, double x_1, double y_1, double constant, vector<pair<double, double>> domain) {
	pair<double, double> pos_temp;

	double temp_max = 0;
	double x_start, x_end, y_start, y_end;
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

		if (abs(x_start - x_end) <= DBL_EPSILON) {
			y_max_point = -(xy*x_start + y_1) / (2 * y_2);
			x_max_point = x_start;

			temp_max = x_2 * x_start*x_start + xy * x_start*y_max_point + y_2 * y_max_point*y_max_point + x_1 * x_start + y_1 * y_max_point + constant;

			if (global_max < temp_max && y_max_point > min(y_start, y_end) && y_max_point <max(y_start, y_end)) {
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
		else {
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

bool check_point(double p_x, double p_y, vector<pair<double, double>>domain) {
	pair<double, double> temp;
	pair<double, double> temp2;
	int check = 0;
	int i = 0;
	for (i = 0; i < domain.size() - 1; i++) {
		temp = pair<double, double>(domain[i + 1].first - domain[i].first, domain[i + 1].second - domain[i].second);
		temp2 = pair<double, double>(p_x - domain[i].first, p_y - domain[i].second);
		if (temp.first*temp2.second - temp2.first*temp.second < -DBL_EPSILON)
			return false;

	}
	temp = pair<double, double>(domain[0].first - domain[i].first, domain[0].second - domain[i].second);
	temp2 = pair<double, double>(p_x - domain[i].first, p_y - domain[i].second);
	if (temp.first*temp2.second - temp2.first*temp.second < -DBL_EPSILON)
		return false;
	return true;
}
