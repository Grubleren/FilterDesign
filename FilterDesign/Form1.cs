using System;
using System.Windows.Forms;
using System.Threading;

namespace JH.Applications
{
    public partial class Form1 : Form
    {
        int fmax;
        int N;

        double om_p;
        double om_s;
        double om_c;
        double a_p;
        double a_s;
        double k;
        double k1;

        public Form1()
        {
            InitializeComponent();

            fmax = 24000;
            N = 2400;

        }

        private void Form1_Load(object sender, EventArgs e)
        {
            Thread thread = new Thread(new ThreadStart(Go));
            thread.Start();
        }

        void Go()
        {
            Complex[] poles;
            Complex[] newpoles;
            Complex[] zeros;
            Complex[] newzeros;

            double eps;
            int n;
            double gain;

            // Specification for LP
            double fp = 10000;
            double fs = 14000;
            a_p = 1;
            a_s = 60;

            double om1 = 2 * Math.PI * fp;
            double om2 = 2 * Math.PI * fs;
            om_p = om1;
            om_s = om2;

            k = om_p / om_s;
            k1 = (Math.Pow(10, a_p / 10) - 1) / (Math.Pow(10, a_s / 10) - 1);

            Butterworth(out n, out poles);

            ButterworthLP(poles, out newpoles, out gain);
            ShowFrf(gain, null, newpoles, 0);

            Chebyshev(out n, out eps, out poles);

            ChebyshevLP(eps, poles, out newpoles, out gain);
            ShowFrf(gain, null, newpoles, 1);

            InvChebyshev(out n, out eps, out zeros, out poles);

            InvChebyshevLP(zeros, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 2);

            // Specification for HP
            fp = 14000;
            fs = 10000;
            a_p = 1;
            a_s = 60;

            om1 = 2 * Math.PI * fp;
            om2 = 2 * Math.PI * fs;
            om_p = 1 / om1;
            om_s = 1 / om2;

            k = om_p / om_s;
            k1 = (Math.Pow(10, a_p / 10) - 1) / (Math.Pow(10, a_s / 10) - 1);

            Butterworth(out n, out poles);

            ButterworthHP(out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 3);

            Chebyshev(out n, out eps, out poles);

            ChebyshevHP(eps, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 4);

            InvChebyshev(out n, out eps, out zeros, out poles);

            InvChebyshevHP(zeros, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 5);

            // Specification for BP
            double fp1 = 10000;
            double fp2 = 14000;
            double fs2 = 16000;
            double fs1 = fp1 * fp2 / fs2;
            a_p = 1;
            a_s = 60;

            om_c = 2 * Math.PI * Math.Sqrt(fp1 * fp2);
            om1 = 2 * Math.PI * fp2;
            om2 = 2 * Math.PI * fs2;
            om_p = (om1 / om_c - om_c / om1);
            om_s = (om2 / om_c - om_c / om2);

            k = om_p / om_s;
            k1 = (Math.Pow(10, a_p / 10) - 1) / (Math.Pow(10, a_s / 10) - 1);

            Butterworth(out n, out poles);

            ButterworthBP(out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 6);

            Chebyshev(out n, out eps, out poles);

            ChebyshevBP(eps, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 7);

            InvChebyshev(out n, out eps, out zeros, out poles);

            InvChebyshevBP(zeros, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 8);

            // Specification for digital BP
            double f_s = 48000;
            fp1 = f_s/Math.PI*Math.Tan(Math.PI*10000/f_s);
            fp2 = f_s / Math.PI * Math.Tan(Math.PI * 14000 / f_s);
            fs2 = f_s / Math.PI * Math.Tan(Math.PI * 16000 / f_s);
            fs1 = fp1 * fp2 / fs2;
            a_p = 1;
            a_s = 60;

            om_c = 2 * Math.PI * Math.Sqrt(fp1 * fp2);
            om1 = 2 * Math.PI * fp2;
            om2 = 2 * Math.PI * fs2;
            om_p = (om1 / om_c - om_c / om1);
            om_s = (om2 / om_c - om_c / om2);

            k = om_p / om_s;
            k1 = (Math.Pow(10, a_p / 10) - 1) / (Math.Pow(10, a_s / 10) - 1);

            InvChebyshev(out n, out eps, out zeros, out poles);

            InvChebyshevBP(zeros, out newzeros, poles, out newpoles, out gain);
            ShowFrf(gain, newzeros, newpoles, 9);
            Complex[] dzeros;
            Complex[] dpoles;
            BilinearTransform(f_s, newzeros, out dzeros, newpoles, out dpoles);
            gain = DGain(1.0, new Complex(Math.Cos(Math.PI / 2), Math.Sin(Math.PI / 2)), dzeros, dpoles);
            ShowDFrf(gain, dzeros, dpoles, 9);

        }

        void Butterworth(out int n, out Complex[] poles)
        {
            n = (int)Math.Ceiling(0.5 * Math.Log10(k1) / Math.Log10(k));
            poles = new Complex[n];

            for (int i = 0; i < poles.Length; i++)
                poles[i] = new Complex(-Math.Sin(Math.PI * (2 * i + 1) / (2 * n)), Math.Cos(Math.PI * (2 * i + 1) / (2 * n)));
        }

        void ButterworthLP(Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            Complex[] zeros = null;
            Complex[] newzeros = null;

            double om0 = om_p * Math.Pow(Math.Pow(10, a_p / 10) - 1, -1.0 / (2 * poles.Length));

            PoleZeroLP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0, 0.0, newzeros, newpoles);
        }

        void ButterworthHP(out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            Complex[] zeros = null;
            newzeros = new Complex[poles.Length];

            double om0 = om_p * Math.Pow(Math.Pow(10, a_p / 10) - 1, -1.0 / (2 * poles.Length));

            PoleZeroHP(om0, zeros, newzeros, poles, newpoles);

            gain = 1.0;
        }

        void ButterworthBP(out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[2 * poles.Length];
            Complex[] zeros = null;
            newzeros = new Complex[poles.Length];
            
            double om0 = om_p * Math.Pow(Math.Pow(10, a_p / 10) - 1, -1.0 / (2 * poles.Length));

            PoleZeroBP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0, om_c, newzeros, newpoles);
        }

        void Chebyshev(out int n, out double eps, out Complex[] poles)
        {
            n = (int)Math.Ceiling(Math.Log10(2 / Math.Sqrt(k1)) / Math.Log10(1 / k + Math.Sqrt(Math.Pow(k, -2) - 1)));
            eps = Math.Sqrt(Math.Pow(10, a_p / 10) - 1);
            poles = new Complex[n];

            double u = 1.0 / n * ArcSinh(1 / eps);
            for (int i = 0; i < poles.Length; i++)
                poles[i] = new Complex(-Math.Sinh(u) * Math.Sin(Math.PI * (2 * i + 1) / (2 * n)), Math.Cosh(u) * Math.Cos(Math.PI * (2 * i + 1) / (2 * n)));
        }

        void ChebyshevLP(double eps, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            Complex[] zeros = null;
            Complex[] newzeros = null;

            double om0 = om_p;

            PoleZeroLP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0 / Math.Sqrt(1 + eps * eps), 0.0, newzeros, newpoles);

        }

        void ChebyshevHP(double eps, out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            Complex[] zeros = null;
            newzeros = new Complex[poles.Length];

            double om0 = om_p;

            PoleZeroHP(om0, zeros, newzeros, poles, newpoles);

            gain = 1.0 / Math.Sqrt(1 + eps * eps);
        }

        void ChebyshevBP(double eps, out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[2 * poles.Length];
            Complex[] zeros = null;
            newzeros = new Complex[poles.Length];
            
            double om0 = om_p;

            PoleZeroBP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0, om_c, newzeros, newpoles);
        }

        void InvChebyshev(out int n, out double eps, out Complex[] zeros, out Complex[] poles)
        {
            n = (int)Math.Ceiling(Math.Log10(2 / Math.Sqrt(k1)) / Math.Log10(1 / k + Math.Sqrt(Math.Pow(k, -2) - 1)));
            eps = Math.Pow(10, -a_s / 20);

            poles = new Complex[n];

            if ((n & 1) == 1)
            {
                zeros = new Complex[n - 1];

                int j = 0;
                for (int i = 0; i < zeros.Length; i++)
                {
                    zeros[i] = 1 / new Complex(0, Math.Cos(Math.PI * (2 * j + 1) / (2 * n)));
                    if (i == n / 2 - 1)
                        j++;
                    j++;
                }
            }
            else
            {
                zeros = new Complex[n];

                for (int i = 0; i < zeros.Length; i++)
                    zeros[i] = 1 / new Complex(0, Math.Cos(Math.PI * (2 * i + 1) / (2 * n)));
            }

            double u = 1.0 / n * ArcSinh(1 / eps);
            for (int i = 0; i < poles.Length; i++)
            {
                poles[i] = 1 / new Complex(-Math.Sinh(u) * Math.Sin(Math.PI * (2 * i + 1) / (2 * n)), Math.Cosh(u) * Math.Cos(Math.PI * (2 * i + 1) / (2 * n)));
            }
        }

        void InvChebyshevLP(Complex[] zeros, out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            newzeros = new Complex[zeros.Length];

            double om0 = om_s;

            PoleZeroLP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0, 0.0, newzeros, newpoles);
        }

        void InvChebyshevHP(Complex[] zeros, out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[poles.Length];
            newzeros = new Complex[zeros.Length + poles.Length - zeros.Length];

            double om0 = om_s;

            PoleZeroHP(om0, zeros, newzeros, poles, newpoles);

            gain = 1.0;
        }

        void InvChebyshevBP(Complex[] zeros, out Complex[] newzeros, Complex[] poles, out Complex[] newpoles, out double gain)
        {
            newpoles = new Complex[2 * poles.Length];
            newzeros = new Complex[2 * zeros.Length + poles.Length - zeros.Length];
            
            double om0 = om_s;

            PoleZeroBP(om0, zeros, newzeros, poles, newpoles);

            gain = Gain(1.0, om_c, newzeros, newpoles);
        }

        void PoleZeroLP(double om0, Complex[] zeros, Complex[] newzeros, Complex[] poles, Complex[] newpoles)
        {
            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                    newzeros[i] = om0 * zeros[i];
            }

            for (int i = 0; i < poles.Length; i++)
                newpoles[i] = om0 * poles[i];
        }

        void PoleZeroHP(double om0, Complex[] zeros, Complex[] newzeros, Complex[] poles, Complex[] newpoles)
        {
            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                    newzeros[i] = 1 / (om0 * zeros[i]);
            }

            for (int i = 0; i < poles.Length; i++)
                newpoles[i] = 1 / (om0 * poles[i]);
        }

        void PoleZeroBP(double om0, Complex[] zeros, Complex[] newzeros, Complex[] poles, Complex[] newpoles)
        {
            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                {
                    newzeros[2 * i] = BPSqrt(om0 * zeros[i], 1);
                    newzeros[2 * i + 1] = BPSqrt(om0 * zeros[i], -1);
                }
            }

            for (int i = 0; i < poles.Length; i++)
            {
                newpoles[2 * i] = BPSqrt(om0 * poles[i], 1);
                newpoles[2 * i + 1] = BPSqrt(om0 * poles[i], -1);
            }
        }

        Complex BPSqrt(Complex root, int sign)
        {
            Complex a = root * root - 4;
            double asq = Complex.AbsSqr(a);
            double angle = Math.Atan2(a.im, a.re);
            Complex sqrt = Math.Pow(asq, 0.25) * new Complex(Math.Cos(angle / 2), Math.Sin(angle / 2));
            return om_c * 0.5 * (root + sign * sqrt);
        }

        double ArcSinh(double y)
        {
            double x;
            x = Math.Log(2 * y);
            double delta = double.MaxValue;

            while (Math.Abs(delta) > 0.0000000001)
            {
                delta = (y - Math.Sinh(x)) / Math.Cosh(x);

                x += delta;

            }

            return x;
        }

        double Gain(double gain, double omref, Complex[] zeros, Complex[] poles)
        {
            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                    gain /= Math.Sqrt(Complex.AbsSqr(new Complex(0, omref) - zeros[i]));
            }
            for (int i = 0; i < poles.Length; i++)
                gain *= Math.Sqrt(Complex.AbsSqr(new Complex(0, omref) - poles[i]));

            return gain;
        }

        double DGain(double gain, Complex omref, Complex[] zeros, Complex[] poles)
        {
            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                    gain /= Math.Sqrt(Complex.AbsSqr(omref - zeros[i]));
            }
            for (int i = 0; i < poles.Length; i++)
                gain *= Math.Sqrt(Complex.AbsSqr(omref - poles[i]));

            return gain;
        }

        void BilinearTransform(double fs, Complex[] zeros, out Complex[] dzeros, Complex[] poles, out Complex[] dpoles)
        {
            dzeros = new Complex[zeros.Length];
            dpoles = new Complex[poles.Length];

            if (zeros != null)
            {
                for (int i = 0; i < zeros.Length; i++)
                    dzeros[i] = (1 + zeros[i] / 2 / fs) / (1 - zeros[i] / 2 / fs);
            }

            for (int i = 0; i < poles.Length; i++)
                dpoles[i] = (1 + poles[i] / 2 / fs) / (1 - poles[i] / 2 / fs);

        }

        Complex[] FRF(double gain, Complex[] zeros, Complex[] poles)
        {
            Complex[] H = new Complex[N];

            if (zeros != null)
            {
                for (int i = 0; i < N; i++)
                {
                    Complex jom = new Complex(0, 2 * Math.PI * i / N * fmax);
                    H[i] = new Complex(gain, 0);
                    for (int j = 0; j < zeros.Length; j++)
                    {
                        H[i] *= (jom - zeros[j]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < N; i++)
                {
                    H[i] = new Complex(gain, 0);
                }
            }

            for (int i = 0; i < N; i++)
            {
                Complex jom = new Complex(0, 2 * Math.PI * i / N * fmax);
                for (int j = 0; j < poles.Length; j++)
                {
                    H[i] *= 1 / (jom - poles[j]);
                }
            }

            return H;
        }

        Complex[] DFRF(double gain, Complex[] zeros, Complex[] poles)
        {
            Complex[] H = new Complex[N];

            if (zeros != null)
            {
                for (int i = 0; i < N; i++)
                {
                    Complex jom = new Complex(Math.Cos(Math.PI * i / N), Math.Sin(Math.PI * i / N));
                    H[i] = new Complex(gain, 0);
                    for (int j = 0; j < zeros.Length; j++)
                    {
                        H[i] *= (jom - zeros[j]);
                    }
                }
            }
            else
            {
                for (int i = 0; i < N; i++)
                {
                    H[i] = new Complex(gain, 0);
                }
            }

            for (int i = 0; i < N; i++)
            {
                Complex jom = new Complex(Math.Cos(Math.PI * i / N), Math.Sin(Math.PI * i / N));
                for (int j = 0; j < poles.Length; j++)
                {
                    H[i] *= 1 / (jom - poles[j]);
                }
            }

            return H;
        }

        float[] Power(Complex[] H)
        {
            float[] power = new float[N];

            for (int i = 0; i < N; i++)
            {
                power[i] = (float)Complex.AbsSqr(H[i]);
            }

            return power;
        }

        void ShowFrf(double gain, Complex[] zeros, Complex[] poles, int nGraph)
        {
            Complex[] H = FRF(gain, zeros, poles);
            float[] P = Power(H);

        }

        void ShowDFrf(double gain, Complex[] zeros, Complex[] poles, int nGraph)
        {
            Complex[] H = DFRF(gain, zeros, poles);
            float[] P = Power(H);

        }

    }
}