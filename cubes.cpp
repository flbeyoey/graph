#include <bits/stdc++.h>
using namespace std;

enum RGBColer
{
	Red = 0,
	Green,
	Bule
};

// #define X dot[0]
// #define Y dot[1]
// #define Z dot[2]
// #define W dot[3]
const int kTriangle_Numbers[15][3] = {{0, 1, 4}, {5, 1, 4}, {6, 4, 7}, {5, 4, 7}, {6, 2, 4}, {0, 2, 4}, {7, 3, 6}, {2, 3, 6}, {5, 1, 7}, {3, 1, 7}, {0, 1, 2}, {0, 2, 3}};
const int kTriangle_color[15][3] = {{200, 0, 0}, {0, 200, 0}, {0, 0, 200}, {200, 200, 0}, {200, 0, 200}, {0, 200, 200}, {100, 200, 0}, {200, 100, 0}, {100, 0, 200}, {200, 0, 100}, {0, 200, 100}, {0, 100, 200}};
const int kScreenLength = 1e3;
const double kesp = 0.000001;

int ScreenBuffer[kScreenLength + 1][kScreenLength + 1][3];
double Mindeepth[kScreenLength + 1][kScreenLength + 1], MXdeepth = -1, MNdeepth = 999;
int TriangleForPixel[kScreenLength + 1][kScreenLength + 1];
int CubeNumber;

struct Matrix4
{
	double element[4][4];
	// 4 X() 4 times 4 X() 4
	Matrix4 operator^(const Matrix4 &matrix)
	{
		Matrix4 ret = {0};
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int k = 0; k < 4; k++)
					ret.element[i][j] += element[i][k] * matrix.element[k][j];
		return ret;
	}
};

struct Node
{
	double dot[4];

	double &X()
	{
		return dot[0];
	}
	double &Y()
	{
		return dot[1];
	}
	double &Z()
	{
		return dot[2];
	}
	double &W()
	{
		return dot[3];
	}
	double X() const
	{
		return dot[0];
	}
	double Y() const
	{
		return dot[1];
	}
	double Z() const
	{
		return dot[2];
	}
	double W() const
	{
		return dot[3];
	}

	void Unitization()
	{
		double len = sqrt(X() * X() + Y() * Y() + Z() * Z());
		(*this) = (*this) / len;
	}
	void Normalized()
	{
		if (!dot[3])
			return;
		X() /= W();
		Y() /= W();
		Z() /= W();
		W() = 1;
	}
	Node operator-(Node node)
	{
		return Node{X() - node.X(), Y() - node.Y(), Z() - node.Z(), W() - node.W()};
	}
	Node operator+(Node node)
	{
		Node point = Node{X() + node.X(), Y() + node.Y(), Z() + node.Z(), W() + node.W()};
		point.Normalized();
		return point;
	}
	// number times
	Node operator*(double Double)
	{
		return Node{X() * Double, Y() * Double, Z() * Double, 0};
	}
	Node operator/(double Double)
	{
		return (*this) * (1 / Double);
	}
	// cross times
	Node operator^(Node node)
	{
		return Node{Y() * node.Z() - Z() * node.Y(),
					Z() * node.X() - X() * node.Z(),
					X() * node.Y() - Y() * node.X(), 0};
	}
	// dot times
	double operator&(Node node)
	{
		return X() * node.X() + Y() * node.Y() + Z() * node.Z();
	}
	Node operator|(Matrix4 matrix)
	{
		Node ret = {0};
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				ret.dot[i] += dot[j] * matrix.element[i][j];
		ret.Normalized();
		return ret;
	}
	// 2D cross times
	double operator%(Node node)
	{
		return X() * node.Y() - Y() * node.X();
	}
	void output()
	{
		printf("%lf %lf %lf %lf\n", X(), Y(), Z(), W());
	}
	double length()
	{
		return sqrt(X() * X() + Y() * Y() + Z() * Z());
	}

	// i suggest that the projection plane is z=2

	Node MVP_P()
	{
		Matrix4 matrix = {0};
		matrix.element[1][1] = matrix.element[2][2] = matrix.element[0][0] = 1;
		matrix.element[3][3] = Z() / 2.0;
		Node ret = (*this) | matrix;
		ret.Z() = 2;
		return ret;
	}
};

struct Triangle
{
	Node a, b, c;
	// This is a 2D check
	bool CheckInTheTriangle(Node node)
	{
		Node veca = node - a, vecb = node - b, vecc = node - c;
		Node ab = b - a, bc = c - b, ca = a - c;
		return ((ab % veca) >= 0 && (bc % vecb) >= 0 && (ca % vecc) >= 0) ||
			   ((ab % veca) <= 0 && (bc % vecb) <= 0 && (ca % vecc) <= 0);
	}
	void TransferIntoProjectionPlane()
	{
		a = a.MVP_P();
		b = b.MVP_P();
		c = c.MVP_P();
	}
};

struct Transform
{
	Node origin;
	Node x_axis, y_axis, z_axis;
	double proportionx, proportiony, proportionz;
};

// some important node
Transform Camera, Cube[105], World;

void MVP_MV(Triangle *triangles,Transform cube)
{
	Node node[8] = {0};
	Node xasia = Node{1, 0, 0, 0}, yasia = Node{0, 1, 0, 0}, zasia = Node{0, 0, 1, 0};
	int tot = 0;

	// M change
	Matrix4 sizeM = {0}, placeM = {0};
	Matrix4 angleM, synthesisM;
	sizeM.element[0][0] = cube.proportionx;
	sizeM.element[1][1] = cube.proportiony;
	sizeM.element[2][2] = cube.proportionz;
	sizeM.element[3][3] = 1;

	placeM.element[0][3] = cube.origin.X();
	placeM.element[1][3] = cube.origin.Y();
	placeM.element[2][3] = cube.origin.Z();
	placeM.element[0][0] = placeM.element[1][1] = placeM.element[2][2] = placeM.element[3][3] = 1;

	angleM = Matrix4{{{cube.x_axis.X(), cube.y_axis.X(), cube.z_axis.X()},
					  {cube.x_axis.Y(), cube.y_axis.Y(), cube.z_axis.Y()},
					  {cube.x_axis.Z(), cube.y_axis.Z(), cube.z_axis.Z()},
					  {0, 0, 0, 1}}};

	synthesisM = placeM ^ (sizeM ^ angleM);

	for (int i = 0; i <= 1; i++)
		for (int j = 0; j <= 1; j++)
			for (int k = 0; k <= 1; k++)
			{
				int a = sqrt(45);
				node[tot] = Node{{0, 0, 0, 1}} + (xasia * i) + (yasia * j) + (zasia * k);
				node[tot] = node[tot] | synthesisM;
				// printf("%lf %lf %lf\n", node[tot].X(), node[tot].Y(), node[tot].Z());
				tot++;
			}
	tot = 0;
	// V change
	sizeM.element[0][0] = Camera.proportionx;
	sizeM.element[1][1] = Camera.proportiony;
	sizeM.element[2][2] = Camera.proportionz;
	sizeM.element[3][3] = 1;

	angleM = Matrix4{{{Camera.x_axis.X(), Camera.x_axis.Y(), Camera.x_axis.Z()},
					  {Camera.y_axis.X(), Camera.y_axis.Y(), Camera.y_axis.Z()},
					  {-Camera.z_axis.X(), -Camera.z_axis.Y(), -Camera.z_axis.Z()},
					  {0, 0, 0, 1}}};

	Node newcamera = Camera.origin | angleM;

	placeM.element[0][3] = -newcamera.X();
	placeM.element[1][3] = -newcamera.Y();
	placeM.element[2][3] = -newcamera.Z();
	placeM.element[0][0] = placeM.element[1][1] = placeM.element[2][2] = placeM.element[3][3] = 1;

	synthesisM = placeM ^ (sizeM ^ angleM);
	for (int i = 0; i < 8; i++)
	{
		node[i] = node[i] | synthesisM;
		// printf("%lf %lf %lf\n", node[i].X(), node[i].Y(), node[i].Z());
	}

	for (int i = 0; i < 12; i++)
		triangles[i] = Triangle{node[kTriangle_Numbers[i][0]],
								node[kTriangle_Numbers[i][1]],
								node[kTriangle_Numbers[i][2]]};
}

// at the first version, i just used a simple way to try to finish it
// and i would add more details at the next version

// there are two founction i need to solve
// first i should get a smaller area to make the program more efficient

// Next i would set a fov parameter and set two projection planes to cut some invisible graphics

void PerspectiveTransformation(Triangle *triangles,int colorfix)
{
	for (int u = 0; u < 12; u++)
	{
		//		if(u == 0) continue;
		//		if(u == 1) continue;
		//		if(u == 2) continue;
		//		if(u == 3) continue;
		//		if(u == 4) continue;
		//		if(u == 5) continue;
		//		if(u == 6) continue;
		//		if(u == 7) continue;
		//		if(u == 8) continue;
		//		if(u == 9) continue;
		//		if(u == 10) continue;
		//		if(u == 11) continue;

		Triangle triangle = triangles[u];
		triangle.TransferIntoProjectionPlane();
		double area = fabs((triangle.b - triangle.a) % (triangle.c - triangle.a));

		// printf("%lf %lf %lf\n", triangles[u].a.X(), triangles[u].a.Y(), triangles[u].a.Z());
		// printf("%lf %lf %lf\n", triangles[u].b.X(), triangles[u].b.Y(), triangles[u].b.Z());
		// printf("%lf %lf %lf\n", triangles[u].c.X(), triangles[u].c.Y(), triangles[u].c.Z());

		// printf("%lf %lf %lf\n", triangle.a.X(), triangle.a.Y(), triangle.a.Z());
		// printf("%lf %lf %lf\n", triangle.b.X(), triangle.b.Y(), triangle.b.Z());
		// printf("%lf %lf %lf\n", triangle.c.X(), triangle.c.Y(), triangle.c.Z());

		// puts("");

		for (int i = 1; i <= 1000; i++)
		{
			for (int j = 1; j <= 1000; j++)
			{
				Node projection = Node{(double)(i - 500) / 1000.0, (double)(j - 500) / 1000.0, 2, 1};
				if (triangle.CheckInTheTriangle(projection))
				{
					double gamma = fabs((projection - triangle.a) % (triangle.b - triangle.a)) / (area),
						   alpha = fabs((projection - triangle.b) % (triangle.c - triangle.b)) / (area),
						   beta = 1 - alpha - gamma,
						   deepth = alpha * triangles[u].a.Z() +
									beta * triangles[u].b.Z() +
									gamma * triangles[u].c.Z();

					bool color_change_tag = false;
					if (Mindeepth[i][j] <= kesp || Mindeepth[i][j] > deepth)
					{
						color_change_tag = true;
						Mindeepth[i][j] = deepth;
						TriangleForPixel[i][j] = u+colorfix*(rand()%2);
					}
				}
			}
		}
	}
}

// get the world coordinate of the camera and the lower left corner of a cube
// get the angle of view of the camera
// get the X() axis and Y() axis of the cube
// the X() axis and Y() axis of the camera are random
// get the proportion of the cube's and camera's Coordinate System

void GetData()
{
	Node &camera = Camera.origin;
	Node lookat;
	scanf("%lf%lf%lf", &camera.X(), &camera.Y(), &camera.Z());

	scanf("%lf%lf%lf", &lookat.X(), &lookat.Y(), &lookat.Z());
	scanf("%lf%lf%lf", &Camera.proportionx, &Camera.proportiony, &Camera.proportionz);
	camera.W() = lookat.W() = 1;
	Camera.z_axis = camera - lookat;
	Camera.z_axis.Unitization();
	// the X() axis and Y() axis of the camera are random
	Node tmp = Node{1, 0, 0, 0};
	Camera.x_axis = Camera.z_axis ^ tmp;
	Camera.y_axis = Camera.x_axis ^ Camera.z_axis;
	Camera.x_axis.Unitization();
	Camera.y_axis.Unitization();

	scanf(" %d ", &CubeNumber);
	// get the cube data
	for (int i = 0; i < CubeNumber; i++)
	{
		scanf("%lf%lf%lf", &Cube[i].origin.X(), &Cube[i].origin.Y(), &Cube[i].origin.Z());
		Cube[i].origin.W() = 1;
		scanf("%lf%lf%lf", &Cube[i].x_axis.X(), &Cube[i].x_axis.Y(), &Cube[i].x_axis.Z());
		scanf("%lf%lf%lf", &Cube[i].y_axis.X(), &Cube[i].y_axis.Y(), &Cube[i].y_axis.Z());
		scanf("%lf%lf%lf", &Cube[i].z_axis.X(), &Cube[i].z_axis.Y(), &Cube[i].z_axis.Z());
		scanf("%lf%lf%lf", &Cube[i].proportionx, &Cube[i].proportiony, &Cube[i].proportionz);
		Cube[i].x_axis.Unitization();
		Cube[i].y_axis.Unitization();
		Cube[i].z_axis.Unitization();
	}
}

void Output()
{

	double _k = 255 / (MXdeepth - MNdeepth);

	for (int i = 1; i <= 1000; i++)
		for (int j = 1; j <= 1000; j++)
		{
			if(Mindeepth[i][j] <= kesp){
				puts("255 255 255");
				continue;
			}
			ScreenBuffer[i][j][Green] = kTriangle_color[TriangleForPixel[i][j]][2];
			ScreenBuffer[i][j][Red] = kTriangle_color[TriangleForPixel[i][j]][0];
			ScreenBuffer[i][j][Bule] = kTriangle_color[TriangleForPixel[i][j]][1];
			// assert(0<=ScreenBuffer[i][j][Green] && ScreenBuffer[i][j][Green] <=255);
			printf("%d %d %d\n", ScreenBuffer[i][j][Red],
				   ScreenBuffer[i][j][Green],
				   ScreenBuffer[i][j][Bule]);
		}
}

int main()
{

	freopen("input", "r", stdin);
	freopen("output.ppm", "w", stdout);
	puts("P3 1000 1000 255");

	World = Transform{Node{0, 0, 0, 1}, Node{1, 0, 0, 0}, Node{0, 1, 0, 0}, Node{0, 0, 1, 0}, 1.0, 1.0, 1.0};
	GetData(); // the variable used in this function are all global variable

	Triangle triangles[105][105] = {0};

	for(int i=0;i<CubeNumber;i++){
		MVP_MV(triangles[i],Cube[i]);
		PerspectiveTransformation(triangles[i],i);
	} 
	Output();
	return 0;
}
