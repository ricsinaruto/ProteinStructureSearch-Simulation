#pragma once

namespace WindowsFormApplication {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;
	

	/// <summary>
	/// Summary for Form1
	/// </summary>
	public ref class Form1 : public System::Windows::Forms::Form
	{
	public:
		Form1(void)
		{
			InitializeComponent();
			//
			//TODO: Add the constructor code here
			//
		}

	protected:
		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		~Form1()
		{
			if (components)
			{
				delete components;
			}
		}
	private: System::Windows::Forms::MenuStrip^  menuStrip1;
	protected:
	private: System::Windows::Forms::ToolStripMenuItem^  toolStripMenuItem1;
	private: System::Windows::Forms::ToolStripMenuItem^  toolStripMenuItem2;



	private: System::Windows::Forms::GroupBox^  groupBox1;
	private: System::Windows::Forms::TextBox^  textBox8;
	private: System::Windows::Forms::TextBox^  textBox7;
	private: System::Windows::Forms::Label^  label8;
	private: System::Windows::Forms::Label^  label7;
	private: System::Windows::Forms::TextBox^  textBox6;
	private: System::Windows::Forms::TextBox^  textBox5;
	private: System::Windows::Forms::TextBox^  textBox4;
	private: System::Windows::Forms::Label^  label6;
	private: System::Windows::Forms::Label^  label5;
	private: System::Windows::Forms::Label^  label4;
	private: System::Windows::Forms::TextBox^  textBox3;
	private: System::Windows::Forms::TextBox^  textBox2;
	private: System::Windows::Forms::TextBox^  textBox1;
	private: System::Windows::Forms::Label^  label3;
	private: System::Windows::Forms::Label^  label2;
	private: System::Windows::Forms::Label^  label1;
	private: System::Windows::Forms::Timer^  timer1;
	private: System::Windows::Forms::PictureBox^  pictureBox1;
	private: System::ComponentModel::IContainer^  components;

	private:
		/// <summary>
		/// Required designer variable.
		/// </summary>
		

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent(void)
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->menuStrip1 = (gcnew System::Windows::Forms::MenuStrip());
			this->toolStripMenuItem1 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->toolStripMenuItem2 = (gcnew System::Windows::Forms::ToolStripMenuItem());
			this->groupBox1 = (gcnew System::Windows::Forms::GroupBox());
			this->textBox8 = (gcnew System::Windows::Forms::TextBox());
			this->textBox7 = (gcnew System::Windows::Forms::TextBox());
			this->label8 = (gcnew System::Windows::Forms::Label());
			this->label7 = (gcnew System::Windows::Forms::Label());
			this->textBox6 = (gcnew System::Windows::Forms::TextBox());
			this->textBox5 = (gcnew System::Windows::Forms::TextBox());
			this->textBox4 = (gcnew System::Windows::Forms::TextBox());
			this->label6 = (gcnew System::Windows::Forms::Label());
			this->label5 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->textBox3 = (gcnew System::Windows::Forms::TextBox());
			this->textBox2 = (gcnew System::Windows::Forms::TextBox());
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->label1 = (gcnew System::Windows::Forms::Label());
			this->timer1 = (gcnew System::Windows::Forms::Timer(this->components));
			this->pictureBox1 = (gcnew System::Windows::Forms::PictureBox());
			this->menuStrip1->SuspendLayout();
			this->groupBox1->SuspendLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->BeginInit();
			this->SuspendLayout();
			// 
			// menuStrip1
			// 
			this->menuStrip1->Items->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->toolStripMenuItem1 });
			this->menuStrip1->Location = System::Drawing::Point(0, 0);
			this->menuStrip1->Name = L"menuStrip1";
			this->menuStrip1->Size = System::Drawing::Size(384, 24);
			this->menuStrip1->TabIndex = 0;
			this->menuStrip1->Text = L"menuStrip1";
			// 
			// toolStripMenuItem1
			// 
			this->toolStripMenuItem1->DropDownItems->AddRange(gcnew cli::array< System::Windows::Forms::ToolStripItem^  >(1) { this->toolStripMenuItem2 });
			this->toolStripMenuItem1->Name = L"toolStripMenuItem1";
			this->toolStripMenuItem1->Size = System::Drawing::Size(125, 20);
			this->toolStripMenuItem1->Text = L"toolStripMenuItem1";
			// 
			// toolStripMenuItem2
			// 
			this->toolStripMenuItem2->Name = L"toolStripMenuItem2";
			this->toolStripMenuItem2->Size = System::Drawing::Size(180, 22);
			this->toolStripMenuItem2->Text = L"toolStripMenuItem2";
			this->toolStripMenuItem2->Click += gcnew System::EventHandler(this, &Form1::toolStripMenuItem2_Click);
			// 
			// groupBox1
			// 
			this->groupBox1->Controls->Add(this->textBox8);
			this->groupBox1->Controls->Add(this->textBox7);
			this->groupBox1->Controls->Add(this->label8);
			this->groupBox1->Controls->Add(this->label7);
			this->groupBox1->Controls->Add(this->textBox6);
			this->groupBox1->Controls->Add(this->textBox5);
			this->groupBox1->Controls->Add(this->textBox4);
			this->groupBox1->Controls->Add(this->label6);
			this->groupBox1->Controls->Add(this->label5);
			this->groupBox1->Controls->Add(this->label4);
			this->groupBox1->Controls->Add(this->textBox3);
			this->groupBox1->Controls->Add(this->textBox2);
			this->groupBox1->Controls->Add(this->textBox1);
			this->groupBox1->Controls->Add(this->label3);
			this->groupBox1->Controls->Add(this->label2);
			this->groupBox1->Controls->Add(this->label1);
			this->groupBox1->Location = System::Drawing::Point(12, 36);
			this->groupBox1->Name = L"groupBox1";
			this->groupBox1->Size = System::Drawing::Size(259, 213);
			this->groupBox1->TabIndex = 1;
			this->groupBox1->TabStop = false;
			this->groupBox1->Text = L"groupBox1";
			// 
			// textBox8
			// 
			this->textBox8->Location = System::Drawing::Point(82, 175);
			this->textBox8->Name = L"textBox8";
			this->textBox8->Size = System::Drawing::Size(63, 20);
			this->textBox8->TabIndex = 16;
			// 
			// textBox7
			// 
			this->textBox7->Location = System::Drawing::Point(13, 176);
			this->textBox7->Name = L"textBox7";
			this->textBox7->Size = System::Drawing::Size(57, 20);
			this->textBox7->TabIndex = 15;
			// 
			// label8
			// 
			this->label8->AutoSize = true;
			this->label8->Location = System::Drawing::Point(79, 159);
			this->label8->Name = L"label8";
			this->label8->Size = System::Drawing::Size(35, 13);
			this->label8->TabIndex = 13;
			this->label8->Text = L"label8";
			// 
			// label7
			// 
			this->label7->AutoSize = true;
			this->label7->Location = System::Drawing::Point(10, 159);
			this->label7->Name = L"label7";
			this->label7->Size = System::Drawing::Size(35, 13);
			this->label7->TabIndex = 12;
			this->label7->Text = L"label7";
			// 
			// textBox6
			// 
			this->textBox6->Location = System::Drawing::Point(167, 132);
			this->textBox6->Name = L"textBox6";
			this->textBox6->Size = System::Drawing::Size(55, 20);
			this->textBox6->TabIndex = 11;
			// 
			// textBox5
			// 
			this->textBox5->Location = System::Drawing::Point(79, 132);
			this->textBox5->Name = L"textBox5";
			this->textBox5->Size = System::Drawing::Size(66, 20);
			this->textBox5->TabIndex = 10;
			// 
			// textBox4
			// 
			this->textBox4->Location = System::Drawing::Point(13, 132);
			this->textBox4->Name = L"textBox4";
			this->textBox4->Size = System::Drawing::Size(57, 20);
			this->textBox4->TabIndex = 9;
			// 
			// label6
			// 
			this->label6->AutoSize = true;
			this->label6->Location = System::Drawing::Point(164, 101);
			this->label6->Name = L"label6";
			this->label6->Size = System::Drawing::Size(35, 13);
			this->label6->TabIndex = 8;
			this->label6->Text = L"label6";
			// 
			// label5
			// 
			this->label5->AutoSize = true;
			this->label5->Location = System::Drawing::Point(76, 101);
			this->label5->Name = L"label5";
			this->label5->Size = System::Drawing::Size(35, 13);
			this->label5->TabIndex = 7;
			this->label5->Text = L"label5";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(10, 101);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(35, 13);
			this->label4->TabIndex = 6;
			this->label4->Text = L"label4";
			// 
			// textBox3
			// 
			this->textBox3->Location = System::Drawing::Point(164, 60);
			this->textBox3->Name = L"textBox3";
			this->textBox3->Size = System::Drawing::Size(58, 20);
			this->textBox3->TabIndex = 5;
			// 
			// textBox2
			// 
			this->textBox2->Location = System::Drawing::Point(76, 60);
			this->textBox2->Name = L"textBox2";
			this->textBox2->Size = System::Drawing::Size(69, 20);
			this->textBox2->TabIndex = 4;
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(10, 61);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(60, 20);
			this->textBox1->TabIndex = 3;
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(161, 28);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(35, 13);
			this->label3->TabIndex = 2;
			this->label3->Text = L"label3";
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(73, 28);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(35, 13);
			this->label2->TabIndex = 1;
			this->label2->Text = L"label2";
			// 
			// label1
			// 
			this->label1->AutoSize = true;
			this->label1->Location = System::Drawing::Point(7, 29);
			this->label1->Name = L"label1";
			this->label1->Size = System::Drawing::Size(35, 13);
			this->label1->TabIndex = 0;
			this->label1->Text = L"label1";
			// 
			// timer1
			// 
			this->timer1->Tick += gcnew System::EventHandler(this, &Form1::timer1_Tick);
			// 
			// pictureBox1
			// 
			this->pictureBox1->Location = System::Drawing::Point(322, 12);
			this->pictureBox1->Name = L"pictureBox1";
			this->pictureBox1->Size = System::Drawing::Size(50, 50);
			this->pictureBox1->TabIndex = 2;
			this->pictureBox1->TabStop = false;
			this->pictureBox1->MouseDown += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::pictureBox1_MouseDown);
			this->pictureBox1->MouseUp += gcnew System::Windows::Forms::MouseEventHandler(this, &Form1::pictureBox1_MouseUp);
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->ClientSize = System::Drawing::Size(384, 261);
			this->Controls->Add(this->pictureBox1);
			this->Controls->Add(this->groupBox1);
			this->Controls->Add(this->menuStrip1);
			this->MainMenuStrip = this->menuStrip1;
			this->Name = L"Form1";
			this->Text = L"Form1";
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->SizeChanged += gcnew System::EventHandler(this, &Form1::Form1_SizeChanged);
			this->menuStrip1->ResumeLayout(false);
			this->menuStrip1->PerformLayout();
			this->groupBox1->ResumeLayout(false);
			this->groupBox1->PerformLayout();
			(cli::safe_cast<System::ComponentModel::ISupportInitialize^>(this->pictureBox1))->EndInit();
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
		Graphics ^fg;
		Pen^pen;
		int x0, y0, dx, dy;
		double n, qx, qy, qz, alfa, beta, gamma, d;
		bool axon;

		ref struct p2d
		{
			double x, y;
			p2d(double x, double y)
			{
				this->x = x; this->y = y;
			}
		};

		ref struct p3d
		{
			double x, y, z;
			p3d(double x, double y, double z)
			{
				this->x = x; this->y = y; this->z = z;
			}
			double f2r(double fok)
			{
				return fok / 180 * Math::PI;
			}
			//axonometric projection
			p2d^axon(double qx, double qy, double qz, double alfa, double beta, double gamma)
			{
				p2d^p = gcnew p2d(-x*qx*Math::Cos(f2r(alfa)) + y*qy*Math::Cos(f2r(beta)) - z*qz*Math::Sin(f2r(gamma)), //xi
					-x*qx*Math::Sin(f2r(alfa)) - y*qy*Math::Sin(f2r(beta)) + z*qz*Math::Cos(f2r(gamma))); //eta
				return p;
			}
		};
		//point value class
		Point wvport(p2d^vp) //magnifying input data
		{
			Point p;
			p.X = Convert::ToInt32(vp->x*n) + this->ClientRectangle.Width / 2;
			p.Y = this->ClientRectangle.Height / 2 - Convert::ToInt32(vp->y*n);
			return p;
		}
		void line3d(p3d^p13d, p3d^p23d)
		{
			p2d^ p12d, ^p22d;
			Point p1k, p2k;
			
			p12d = p13d->axon(qx, qy, qz, alfa, beta, gamma);
			p22d = p23d->axon(qx, qy, qz, alfa, beta, gamma);

			p1k = wvport(p12d);
			p2k = wvport(p22d);
			fg->DrawLine(pen, p1k, p2k);
		}

		double f(int i, double x,double y, double z, double rx, double ry, double rz)
		{
			//i show which coordinate it is (0,1,2)
			//x,y,z are the coordinates and rx, ry, rz are the rotation deegres
			double m[3][3] = {};
			//Matrix
			m[0][0] = Math::Cos(ry)*Math::Cos(rz) + Math::Sin(rx)*Math::Sin(ry)*Math::Sin(rz);
			m[0][1] = -Math::Cos(ry)*Math::Sin(rz) + Math::Sin(rx)*Math::Sin(ry)*Math::Cos(rz);
			m[0][2] = Math::Cos(rx)*Math::Sin(ry);
			m[1][0] = Math::Cos(rx)*Math::Sin(rz);
			m[1][1] = Math::Cos(rx)*Math::Cos(rz);
			m[1][2] = -Math::Sin(rx);
			m[2][0] = Math::Sin(ry)*-Math::Cos(rz) + Math::Sin(rx)*Math::Cos(ry)*Math::Sin(rz);
			m[2][1] = Math::Sin(ry)*Math::Sin(rz) + Math::Sin(rx)*Math::Cos(ry)*Math::Cos(rz);
			m[2][2] = Math::Cos(rx)*Math::Cos(ry);

			double k;
			//k if the value of the rotation vector at i
			if (i == 0)
			{
				k = m[0][0]*x+m[0][1]*y+m[0][2]*z;
			}
			if (i == 1)
			{
				k = m[1][0] * x + m[1][1] * y + m[1][2] * z;
			}
			if (i == 2)
			{
				k = m[2][0] * x + m[2][1] * y + m[2][2] * z;
			}
			return k;
		}
		
		double kocka(double rx, double ry, double rz,int i)
		{
			p3d^p1, ^p2, ^p3, ^p4, ^p5, ^p6, ^p7, ^p8;

			//parameters from textboxes
			qx = Convert::ToDouble(textBox1->Text);
			qy = Convert::ToDouble(textBox2->Text);
			qz = Convert::ToDouble(textBox3->Text);
			alfa = Convert::ToDouble(textBox4->Text);
			beta = Convert::ToDouble(textBox5->Text);
			gamma = Convert::ToDouble(textBox6->Text);
			n = Convert::ToDouble(textBox7->Text);
			d = Convert::ToDouble(textBox8->Text);

			fg = this->CreateGraphics();
			pen = gcnew Pen(System::Drawing::Color::Blue);
			groupBox1->Visible = false;
			//fg->Clear(System::Drawing::SystemColors::Control);
			// coordinate system
			p1 = gcnew p3d(0, 0, 0);
			p2 = gcnew p3d(2, 0, 0);
			line3d(p1, p2);
			p2 = gcnew p3d(0, 2, 0);
			line3d(p1, p2);
			p2 = gcnew p3d(0, 0, 2);
			line3d(p1, p2);
			if (i==0) pen->Color = System::Drawing::Color::Red;
			else pen->Color = System::Drawing::SystemColors::Control;
			

			//cube points
			p1 = gcnew p3d(f(0, 0, 0, 0, rx, ry, rz), f(1, 0, 0, 0, rx, ry, rz), f(2, 0, 0, 0, rx, ry, rz));
			p2 = gcnew p3d(f(0, 1, 0, 0, rx, ry, rz), f(1, 1, 0, 0, rx, ry, rz), f(2, 1, 0, 0, rx, ry, rz));
			p3 = gcnew p3d(f(0, 0, 1, 0, rx, ry, rz), f(1, 0, 1, 0, rx, ry, rz), f(2, 0, 1, 0, rx, ry, rz));
			p4 = gcnew p3d(f(0, 1, 1, 0, rx, ry, rz), f(1, 1, 1, 0, rx, ry, rz), f(2, 1, 1, 0, rx, ry, rz));
			p5 = gcnew p3d(f(0, 0, 0, 1, rx, ry, rz), f(1, 0, 0, 1, rx, ry, rz), f(2, 0, 0, 1, rx, ry, rz));
			p6 = gcnew p3d(f(0, 1, 0, 1, rx, ry, rz), f(1, 1, 0, 1, rx, ry, rz), f(2, 1, 0, 1, rx, ry, rz));
			p7 = gcnew p3d(f(0, 0, 1, 1, rx, ry, rz), f(1, 0, 1, 1, rx, ry, rz), f(2, 0, 1, 1, rx, ry, rz));
			p8 = gcnew p3d(f(0, 1, 1, 1, rx, ry, rz), f(1, 1, 1, 1, rx, ry, rz), f(2, 1, 1, 1, rx, ry, rz));

			//edges
			line3d(p1, p2);
			line3d(p1, p3);
			line3d(p2, p4);
			line3d(p3, p4);
			line3d(p5, p6);
			line3d(p5, p7);
			line3d(p6, p8);
			line3d(p7, p8);
			line3d(p1, p5);
			line3d(p2, p6);
			line3d(p3, p7);
			line3d(p4, p8);

			return 0;
		}

	private: System::Void Form1_Load(System::Object^  sender, System::EventArgs^  e)
	{
		this->Text = "3D";
		groupBox1->Text = "parameters";
		label1->Text = "qx";
		label2->Text = "qy";
		label3->Text = "qz";
		label4->Text = "alfa";
		label5->Text = "beta";
		label6->Text = "gamma";
		label7->Text = "magnify";
		label8->Text = "central d";
		// cavalier
		textBox1->Text = Convert::ToString(0.5);
		textBox2->Text = "1";
		textBox3->Text = "1";
		textBox4->Text = "45";
		textBox5->Text = "0";
		textBox6->Text = "0";
		// cheating
		textBox7->Text = "300";
		textBox8->Text = "5";
		toolStripMenuItem1->Text = "draw";
		toolStripMenuItem2->Text = "axon";
		
	}
			

	private: System::Void toolStripMenuItem2_Click(System::Object^  sender, System::EventArgs^  e)
	{
		axon = true;
		
	}
	
	
private: System::Void Form1_SizeChanged(System::Object^  sender, System::EventArgs^  e)
{
	pictureBox1->Location = Point(Form1::Width-78, 12);
}
private: System::Void pictureBox1_MouseDown(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	pictureBox1->BackColor = SystemColors::Desktop;
	Cursor->Hide();
	x0 = MousePosition.X;
	y0 = MousePosition.Y;
	timer1->Enabled = true;
	timer1->Interval = 15;
	timer1->Start();
}
private: System::Void pictureBox1_MouseUp(System::Object^  sender, System::Windows::Forms::MouseEventArgs^  e)
{
	pictureBox1->BackColor = SystemColors::Control;
	Cursor->Show();
	timer1->Enabled = false;
}

		 double lolx, loly;
private: System::Void timer1_Tick(System::Object^  sender, System::EventArgs^  e)
{
	dx = MousePosition.X - x0;
	dy = MousePosition.Y - y0;

	fg = this->CreateGraphics();
	fg->Clear(System::Drawing::SystemColors::Control);
	lolx = lolx + dx;
	loly = loly + dy;
	kocka(-lolx/500,loly/500,0,0);
	System::Windows::Forms::Cursor::Position = Point(x0, y0);
}
};
}
