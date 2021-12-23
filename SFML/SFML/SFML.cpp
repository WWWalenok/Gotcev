#include <iostream>
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <fstream>

#include "Math.h"


static std::chrono::steady_clock::time_point start_time = std::chrono::high_resolution_clock::now();
const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

double gausrand(double S = 1, double U = 0)
{
	return sqrt(-2 * log((1 + rand()) / float(RAND_MAX + 1))) * cos(2 * PI * (rand() / float(RAND_MAX))) * S + U;
}

struct PoliLine : public sf::Drawable
{
public:
	sf::Vertex* m_vertices = 0;
	uint32_t count = 0;

	PoliLine(int _count = 0, sf::Vertex* _m_vertices = 0, sf::Color _color = sf::Color(0, 0, 0, 0))
	{
		count = _count;
		if (count <= 0)
		{
			return;
		}
		if (_m_vertices == 0)
		{
			m_vertices = new sf::Vertex[count];
		}
		else
		{
			m_vertices = _m_vertices;
		}

		for (int i = 0; i < count; i++)
		{
			m_vertices[i].color = _color;
		}
	}


	sf::Vertex& operator[](uint32_t i)
	{
		return m_vertices[i % count];
	}

	void Set(sf::Vertex* n_vertices, uint32_t count)
	{
		m_vertices = n_vertices;
		PoliLine::count = count;
	}

	void SetColor(sf::Color color)
	{
		for (int i = 0; i < count; i++)
		{
			m_vertices[i].color = color;
		}
	}

private:
	virtual void draw(sf::RenderTarget& target, sf::RenderStates states) const
	{
		target.draw(m_vertices, count, sf::PrimitiveType::LineStrip, states);
	}
};

struct Plot : public sf::Drawable
{
	sf::Sprite sprite;
	sf::Texture texture;
	sf::Uint8* pixels;

	int H, W;
	Plot(int H, int W)
	{
		this->H = H;
		this->W = W;

		pixels = new sf::Uint8[H * W * 4];
		texture.create(H, W);
	}

#define SetPixel(plot, x, y, r, g, b, a) {sf::Uint8* p = &plot.pixels[((x % H) * W + y % W) * 4];p[0] = r; p[1] = g; p[2] = b; p[3] = a;}

	inline void Ubdate()
	{
		texture.update(pixels);
		sprite.setTexture(texture);
	}

private:


	virtual void draw(sf::RenderTarget& target, sf::RenderStates states) const
	{
		target.draw(sprite);
	}
};



void main()
{
	uint16_t H = 800, W = H;

	std::srand(std::time(0));

	sf::ContextSettings context_setting(0, 0, 0);
	sf::RenderWindow window(sf::VideoMode(H, W), "SFML window", sf::Style::Default, context_setting);
	sf::CircleShape shape(100.f);
	shape.setFillColor(sf::Color::Green);

	int Count = 50;

	PoliLine
		Base(Count, 0, sf::Color(255, 0, 0, 128)),
		BaseExp(Count, 0, sf::Color::Yellow),
		Ideal(Count, 0, sf::Color(255, 255, 255)),
		Exp1(Count, 0, sf::Color::Red),
		Exp2(Count, 0, sf::Color(255, 0, 255)/*Pink*/),
		Ap1(Count, 0, sf::Color::Blue),
		Ap2(Count, 0, sf::Color::Green);

	float
		* ys = new float[Count],
		* xs = new float[Count];

	float ec, lc, ac, as, ps, hs;

	float K = 5, O = H / 2.0;

	//Generate
	{
		float
			dx = 0.2;
		ec = 1 + rand() / float(RAND_MAX);
		lc = 3 + gausrand();
		ac = 1 + gausrand();
		as = 3 + gausrand();
		hs = 2 + gausrand();
		ps = 0 + rand() / float(RAND_MAX) * 2 * PI;

		if (ps > PI)
		{
			ps -= PI;
			as *= -1;
		}

		float ymin = 1000000, ymax = -1000000;

		for (int i = 0; i < Count; i++)
		{
			xs[i] = i * dx + gausrand(dx * 0.1);
			ys[i] = gausrand(lc * 0.02) + lc * exp(-ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac;
			ymin = fmin(ymin, ys[i]);
			ymax = fmax(ymax, ys[i]);
		}

		K = 0.4 * (H) / fmax(fabs(ymax), fabs(ymin));

		//O = 0.8 * H - K * ymin;

		std::ofstream fout("gens.txt");
		for (int i = 0; i < Count; i++)
			fout << xs[i] << ",";
		fout << std::endl;
		for (int i = 0; i < Count; i++)
			fout << ys[i] << ",";
		fout << std::endl;
		fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;

		fout.close();

		for (int i = 0; i < Count; i++)
			Base[i].position.x = (xs[i] - xs[0]) / (xs[Count - 1] - xs[0]) * W * 0.9 + W * 0.05;

		for (int i = 0; i < Count; i++)
			Base[i].position.y = (O - ys[i] * K);
	}

	{
		std::ofstream fout("exp.txt");
		for (int i = 0; i < Count; i++)
			fout << lc * exp(-ec * xs[i]) * (1) + ac << ",";
		fout << std::endl;
		fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;

		fout.close();

		for (int i = 0; i < Count; i++)
			BaseExp[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1) + ac) * K);
	}

	{
		std::ofstream fout("sin.txt");
		for (int i = 0; i < Count; i++)
			fout << as * sin(ps + (xs[i]) * hs) * lc * exp(-ec * xs[i]) << ",";
		fout << std::endl;
		fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;

		fout.close();
	}

	for (int i = 0; i < Count; i++)
		Ideal[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac) * K);

	// First Aprocsimate

	for (int i = 0; i < Count; i++)
		BaseExp[i].position.x =
		Ideal[i].position.x =
		Exp1[i].position.x =
		Exp2[i].position.x =
		Ap1[i].position.x =
		Ap2[i].position.x =
		Base[i].position.x;

	ac = 0;
	{
		float C = 0;
		for (int i = Count - 1; i > Count * 0.75; i--)
		{
			C++;
			ac += ys[i];
		}
		ac /= C;
	}

	hs = 0;
	ps = 0;
	int imax = 0, imin = 0;
	ec = lc = 0;
	{
		int FF = 3;
		for (int i = 0; i < Count; i++)
		{
			bool pbreak = false;
			if (ys[imax] < ys[i])
				imax = i;
			if (ys[imin] > ys[i])
				imin = i;
			if (ys[imin] <= ys[i])
				if (i != imax)
					if (imin != 0)
					{
						if (i - imin > FF)
							if (pbreak)
								break;
					}
					else
					{
						if (ys[imax] >= ys[i])
							if (i - imax > FF)
								imin = i;
					}
			if (ys[imax] >= ys[i])
				if (i != imin)
					if (imax != 0)
					{
						if (i - imax > FF)
							if (pbreak)
								break;
					}
					else
					{
						if (ys[imin] <= ys[i])
							if (i - imin > FF)
								imax = i;
					}
		}

		float cs = 0.5 * (xs[imin] + xs[imax]);
		hs = -PI / (xs[imin] - xs[imax]);
		ps = PI - hs * cs;
		while (ps < 0)
			ps += 2 * PI;
		while (ps > 2 * PI)
			ps -= 2 * PI;

		float x0 = cs, y0 = 0;
		float x1 = cs + fabs(xs[imin] - xs[imax]), y1 = 0;
		int i0, i1;
		for (int i = 1; i < Count; i++)
		{
			if (y0 == 0)
			{
				if (xs[i] > x0)
				{
					y0 = (ys[i] * (xs[i] - x0) + ys[i - 1] * (x0 - xs[i - 1])) / (xs[i] - xs[i - 1]);
					i0 = i - 1;
				}
			}

			if (y1 == 0)
			{
				if (xs[i] > x1)
				{
					y1 = (ys[i] * (xs[i] - x1) + ys[i - 1] * (x1 - xs[i - 1])) / (xs[i] - xs[i - 1]);
					i1 = i - 1;
				}
			}
		}

		y0 -= ac;
		y1 -= ac;
		if (true)
		{
			int s = std::max(0, i1 - 3);
			for (int i = 0; i < 7; i++)
			{
				y1 += ys[s + i];
			}
			y1 /= 8.0;
			if (y1 < 0.00001)
			{
				y0 += 0.00001 - y1;
				y1 = 0.00001;
			}
		}
		ec = -std::log(y0 / y1) / (x0 - x1);
		lc = y0 * exp(ec * x0);
		ec = 1;
		//ac = 0;
		lc = 0.2;

		as = -(ys[imax] - ys[imin]) / lc;

		if (ps > PI)
		{
			ps -= PI;
			as *= -1;
		}
	}
	{
		std::ofstream fout("A1.txt");

		for (int i = 0; i < Count; i++)
			fout << lc * exp(-ec * xs[i]) * (1) + ac << ",";
		fout << std::endl;
		fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;
		fout.close();

		for (int i = 0; i < Count; i++)
			Exp1[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1) + ac) * K);
	}
	{
		std::ofstream fout("A2.txt");
		for (int i = 0; i < Count; i++)
			fout << lc * exp(-ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac << ",";
		fout << std::endl;
		fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;
		fout.close();

		for (int i = 0; i < Count; i++)
			Ap1[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac) * K);
	}

	float BKrit = 100000000000;

	double A = 0;
	auto t_end = std::chrono::high_resolution_clock::now();
	auto t_start = std::chrono::high_resolution_clock::now();

	float fps = 0, k = 0;
	std::ofstream log("log.txt");
	for (uint64_t t = 0; window.isOpen(); t++)
	{
		t_end = std::chrono::high_resolution_clock::now();
		double dt = std::chrono::duration<double, std::milli>(t_end - t_start).count() / 1000.0;
		t_start = t_end;

		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
		// Frame Math
		for (int i = 0; i < 1; i++)
		{
			float dreaf = 0;
#define errF  powf(fmax(0, abs(t) - dreaf), 4);
			{
				const int CC = 10;

				float
					b_ec[CC],
					b_lc[CC],
					b_ac[CC],
					Ks[CC];

				for (int j = 0; j < CC; j++)
				{
					b_ec[j] = gausrand(4, ec);
					b_ac[j] = gausrand(4, ac);
					b_lc[j] = gausrand(4, lc);
					Ks[j] = 0;
					for (int i = 0; i < Count; i++)
					{
						float t = b_lc[j] * exp(-b_ec[j] * xs[i]) * (1) + b_ac[j];
						t = t - ys[i];
						Ks[j] += errF;
					}
				}

				int B = 0;

				for (int i = 0; i < CC; i++)
				{
					if (Ks[B] > Ks[i])
						B = i;
				}

				for (int i = 0; i < 10000; i++)
				{
					int j = i % CC;
					float
						t_ec = gausrand(1, b_ec[j]) * 0.8 + gausrand(1, b_ec[B]) * 0.2,
						t_lc = gausrand(1, b_lc[j]) * 0.8 + gausrand(1, b_lc[B]) * 0.2,
						t_ac = gausrand(1, b_ac[j]) * 0.8 + gausrand(1, b_ac[B]) * 0.2;

					float K = 0;
					for (int i = 0; i < Count; i++)
					{
						float t = t_lc * exp(-t_ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac;
						t = t - ys[i];
						K += errF;
					}
					if (K < Ks[j])
					{
						b_ec[j] = t_ec;
						b_lc[j] = t_lc;
						b_ac[j] = t_ac;

						Ks[j] = K;

						if (K < Ks[B])
							B = j;
					}
				}
				if (Ks[B] < BKrit)
				{
					ec = b_ec[B];
					lc = b_lc[B];
					//ac = b_ac[B];
					BKrit = Ks[B];
				}
			}

			{
				const int CC = 10;

				float
					b_as[CC],
					b_ps[CC],
					b_hs[CC],
					Ks[CC];

				for (int j = 0; j < CC; j++)
				{
					b_as[j] = gausrand(4, as);
					b_ps[j] = gausrand(4, ps);
					b_hs[j] = gausrand(4, hs);

					while (b_ps[j] < 0)
						b_ps[j] += 2 * PI;
					while (b_ps[j] > 2 * PI)
						b_ps[j] -= 2 * PI;
					if (b_ps[j] > PI)
					{
						b_ps[j] -= PI;
						b_as[j] *= -1;
					}
					Ks[j] = 0;
					for (int i = 0; i < Count; i++)
					{
						float t = lc * exp(-ec * xs[i]) * (1 + b_as[j] * sin(b_ps[j] + (xs[i]) * b_hs[j])) + ac;
						t = t - ys[i];
						Ks[j] += errF;
					}
				}

				int B = 0;

				for (int i = 0; i < CC; i++)
				{
					if (Ks[B] > Ks[i])
						B = i;
				}

				for (int i = 0; i < 10000; i++)
				{
					int j = i % CC;
					float
						t_as = gausrand(1, b_as[j]) * 0.8 + gausrand(1, b_as[B]) * 0.2,
						t_ps = gausrand(1, b_ps[j]) * 0.8 + gausrand(1, b_ps[B]) * 0.2,
						t_hs = gausrand(1, b_hs[j]) * 0.8 + gausrand(1, b_hs[B]) * 0.2;
					while (t_ps < 0)
						t_ps += 2 * PI;
					while (t_ps > 2 * PI)
						t_ps -= 2 * PI;
					if (t_ps > PI)
					{
						t_ps -= PI;
						t_as *= -1;
					}
					float K = 0;
					for (int i = 0; i < Count; i++)
					{
						float t = lc * exp(-ec * xs[i]) * (1 + t_as * sin(t_ps + (xs[i]) * t_hs)) + ac;
						t = t - ys[i];
						K += errF;
					}
					if (K < Ks[j])
					{
						b_as[j] = t_as;
						b_ps[j] = t_ps;
						b_hs[j] = t_hs;

						Ks[j] = K;

						if (K < Ks[B])
							B = j;
					}
				}
				if (Ks[B] < BKrit)
				{
					as = b_as[B];
					ps = b_ps[B];
					hs = b_hs[B];
					BKrit = Ks[B];
				}
			}
		}

		{

			for (int i = 0; i < Count; i++)
				Ap2[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1 + as * sin(ps + (xs[i]) * hs)) + ac) * K);
		}

		{
			for (int i = 0; i < Count; i++)
				Exp2[i].position.y = (O - (lc * exp(-ec * xs[i]) * (1) + ac) * K);
		}



		window.clear();


		//window.draw(Ideal);
		//window.draw(BaseExp);
		//window.draw(Exp1);
		//window.draw(Ap1);
		//window.draw(Exp2);
		window.draw(Ap2);


		window.draw(Base);

		log << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs << std::endl;

		fps += 1 / dt;
		k++;
		if (t % 1 == 0)
		{
			std::cout << fps / k << "    \r";
			fps = fps * 0.05;
			k = k * 0.05;
		}
		window.display();
	}

	std::ofstream fout("A.txt");

	for (int i = 0; i < Count; i++)
		fout << lc * exp(-ec * xs[i]) * 1 + ac << ",";
	fout << std::endl;
	fout << ec << ',' << lc << ',' << ac << ',' << as << ',' << ps << ',' << hs;
	fout.close();


}
