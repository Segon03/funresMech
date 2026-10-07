// ============================================================
// fun/motor_okuyama.cpp
//
// Motor de simulacion del modelo mecanistico de Okuyama (2012, 2026)
// compilado en C++ via Rcpp. Se carga desde R con fun/motor_rcpp.R
// (no hace falta llamar a sourceCpp() a mano).
//
// Contiene TRES versiones del simulador:
//
// (A) parasitism_gamma_lnorm_cpp() + rlnorm0_cpp()
//     COPIA LITERAL de Okuyama (2026), Script S5
//     (papers/Okuyama_Supp_Mat/jen70148-sup-0008-datas5/.../Script S5.cpp).
//     No se toca ni una linea: sirve de REFERENCIA.
//
// (A') parasitism_gamma_lnorm_cpp_corr()
//     Igual a (A) salvo un detalle tecnico de Rcpp (RNGScope anidado)
//     que, con las versiones actuales de Rcpp, hace que (A) reutilice
//     numeros aleatorios cuando s > 0. Ver el bloque (A') abajo.
//
// (B) motor_okuyama_cpp()
//     El motor que usan los ajustes de este proyecto. Hace
//     EXACTAMENTE el mismo proceso estocastico que (A), consumiendo los
//     numeros aleatorios en el mismo orden (con corte_saturacion=false
//     da resultados IDENTICOS a (A') con la misma semilla, y a (A)
//     cuando s = 0; ver R/00b_motor_rcpp_pruebas.R, prueba 1, y el
//     HALLAZGO sobre RNGScope descripto en (A') mas abajo). Agrega:
//
//     1. modo = 0 ("okuyama"): igual a (A).
//          - tiempo de busqueda ~ Gamma(forma=k, tasa=k*a*H^z)
//          - huesped elegido con floor(runif(0,H))
//          - tiempo de manipulacion ~ Lognormal con MEDIA h y DESVIO
//            ESTANDAR s en escala natural (rlnorm0 de Okuyama);
//            s = 0 -> manipulacion fija = h
//          - el ensayo termina cuando t > T
//        modo = 1 ("funresMech"): igual a simulate_trial() de
//          funresMech v1.0.4 (R puro), para poder comparar con el
//          motor anterior del paquete:
//          - huesped elegido con R_unif_index() (lo que usa
//            sample.int(x, 1) internamente en R >= 3.6)
//          - Lognormal con meanlog = log(h) - s^2/2, sdlog = s
//            (s es el desvio en escala LOG, no natural)
//          - el ensayo termina cuando t >= T
//        (Con la misma semilla, modo=1 reproduce bit a bit a
//         funresMech:::simulate_distribution(); prueba 2.)
//
//     2. corte_saturacion (EXACTO, no aproxima nada): el resultado de
//        un ensayo es la cantidad de huespedes DISTINTOS parasitados.
//        Una vez que los H huespedes ya fueron parasitados, seguir
//        simulando no puede cambiar ese numero, asi que se corta el
//        ensayo ahi. Esto elimina la "trampa computacional" de h ~ 0
//        con k ~ 0 (o a*H^z enorme): en esos casos el while(t < T)
//        original daba millones de vueltas sin cambiar el resultado
//        (ver docs/MOTOR_RCPP_cambios_y_validacion.md). Con el corte,
//        cada ensayo hace a lo sumo ~H*log(H) vueltas (problema del
//        coleccionista de figuritas). La DISTRIBUCION del resultado
//        es identica; solo cambia cuantos numeros aleatorios se
//        consumen despues de la saturacion, por eso con la misma
//        semilla los ensayos siguientes difieren individualmente
//        pero no en distribucion (prueba 3).
//
//     3. Pequenas optimizaciones sin efecto sobre el resultado:
//        phi y mu de la lognormal se calculan una vez por llamada (no
//        una vez por evento), y no se crea un NumericVector por cada
//        tiempo de manipulacion.
// ============================================================

#include <Rcpp.h>
#include <R_ext/Random.h>   // R_unif_index (modo funresMech)
#include <cmath>
#include <vector>

using namespace Rcpp;

// ------------------------------------------------------------
// (A) COPIA LITERAL de Okuyama (2026) Script S5 - NO MODIFICAR
// ------------------------------------------------------------

// Custom lognormal generator with specified mean and sd
// [[Rcpp::export]]
NumericVector rlnorm0_cpp(int n, double mean, double sd) {
	RNGScope scope;
    if (sd == 0) {
        return NumericVector(n, mean);
    }
    if (sd < 0 || mean <= 0) {
        stop("Infeasible parameters");
    }

    double phi = sqrt(log(1 + (sd / mean) * (sd / mean)));
    double mu = log(mean) - 0.5 * phi * phi;

    NumericVector result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = exp(R::rnorm(mu, phi));
    }
    return result;
}

// Parasitism simulation for a single H value, returning nsim results
// [[Rcpp::export]]
IntegerVector parasitism_gamma_lnorm_cpp(int H, double a, double h, double z, double k, double s, double T, int nsim = 10000) {
	RNGScope scope;
    IntegerVector parasitized_counts(nsim);

    for (int sim = 0; sim < nsim; ++sim) {
        NumericVector h0(H, 0.0);
        double t = 0.0;

        while (t < T) {
            double rate = k * a * pow(H, z);
            double search = R::rgamma(k, 1.0 / rate);
            t += search;
            if (t > T) break;

            int pid = floor(R::runif(0, H));
            h0[pid] += 1;

            double handling_time = (s == 0) ? h : rlnorm0_cpp(1, h, s)[0];
            t += handling_time;
        }

        int count = 0;
        for (int i = 0; i < H; ++i) {
            if (h0[i] > 0) count++;
        }

        parasitized_counts[sim] = count;

        // Yield control every 1000 simulations to prevent freezing
        if (sim % 1000 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }

    return parasitized_counts;
}

// ------------------------------------------------------------
// (A') Script S5 con UNA sola correccion: la lognormal se genera sin
//      abrir un RNGScope anidado.
//
// HALLAZGO (28-sep-2026, ver R/00b prueba 1 y docs/): al menos desde
// Rcpp 1.0.0 (2018) RNGScope ya NO tiene contador: cada RNGScope llama a
// GetRNGstate() al abrirse, que RELEE .Random.seed. Como rlnorm0_cpp()
// abre su propio RNGScope y se llama DENTRO del bucle de
// parasitism_gamma_lnorm_cpp() (que ya tiene uno abierto y todavia no
// escribio el estado con PutRNGstate()), el generador se "rebobina":
// el tiempo de manipulacion reutiliza numeros aleatorios que ya se
// usaron para el tiempo de busqueda y para elegir el huesped. El
// resultado es una distribucion SESGADA siempre que s > 0 (con s = 0
// no se llama a rlnorm0_cpp y el problema no aparece).
// La unica diferencia de (A') con (A) es esa: mismo codigo, mismo
// orden de consumo de numeros aleatorios, pero sin el RNGScope interno.
// ------------------------------------------------------------

static inline double rlnorm0_interno(double mean, double sd) {
    double phi = sqrt(log(1 + (sd / mean) * (sd / mean)));
    double mu = log(mean) - 0.5 * phi * phi;
    return exp(R::rnorm(mu, phi));
}

// [[Rcpp::export]]
IntegerVector parasitism_gamma_lnorm_cpp_corr(int H, double a, double h, double z, double k, double s, double T, int nsim = 10000) {
	RNGScope scope;
    IntegerVector parasitized_counts(nsim);

    for (int sim = 0; sim < nsim; ++sim) {
        NumericVector h0(H, 0.0);
        double t = 0.0;

        while (t < T) {
            double rate = k * a * pow(H, z);
            double search = R::rgamma(k, 1.0 / rate);
            t += search;
            if (t > T) break;

            int pid = floor(R::runif(0, H));
            h0[pid] += 1;

            double handling_time = (s == 0) ? h : rlnorm0_interno(h, s);   // <- unico cambio
            t += handling_time;
        }

        int count = 0;
        for (int i = 0; i < H; ++i) {
            if (h0[i] > 0) count++;
        }

        parasitized_counts[sim] = count;

        if (sim % 1000 == 0) {
            Rcpp::checkUserInterrupt();
        }
    }

    return parasitized_counts;
}

//' Demostracion minima del efecto del RNGScope anidado: devuelve pares
//' (numero del bucle externo, numero generado dentro de una funcion con
//' su propio RNGScope). Con Rcpp >= 1.0 los pares salen repetidos.
static double runif_con_scope_propio() { RNGScope scope; return R::runif(0, 1); }
// [[Rcpp::export]]
NumericMatrix demo_rngscope_anidado(int n) {
    RNGScope scope;
    NumericMatrix out(n, 2);
    for (int i = 0; i < n; ++i) {
        out(i, 0) = R::runif(0, 1);
        out(i, 1) = runif_con_scope_propio();
    }
    return out;
}

// ------------------------------------------------------------
// (B) Motor de este proyecto
// ------------------------------------------------------------

//' Simula nsim ensayos de parasitismo a densidad H y devuelve, para
//' cada uno, la cantidad de huespedes distintos parasitados.
//' modo: 0 = okuyama (Script S5), 1 = funresMech v1.0.4 (R puro)
//' corte_saturacion: true = cortar el ensayo cuando los H huespedes
//'   ya estan parasitados (exacto; evita la trampa computacional).
// [[Rcpp::export]]
IntegerVector motor_okuyama_cpp(int H, double a, double h, double z,
                                double k, double s, double T,
                                int nsim, int modo = 0,
                                bool corte_saturacion = true) {
    RNGScope scope;
    if (H < 1) stop("H debe ser >= 1");
    if (!(a > 0) || !(k > 0) || !(h >= 0) || !(s >= 0) || !(T > 0))
        stop("parametros no validos (se requiere a>0, k>0, h>=0, s>=0, T>0)");

    IntegerVector out(nsim);
    std::vector<unsigned char> tocado(H);

    // Tasa y escala de la Gamma. Se respeta el orden de las
    // operaciones de cada motor original (el redondeo en punto
    // flotante puede diferir en el ultimo bit segun el orden).
    double escala_gamma;
    if (modo == 0) {
        double rate = k * a * std::pow((double) H, z);        // Script S5
        escala_gamma = 1.0 / rate;
    } else {
        double lambda = a * R_pow((double) H, z);             // simulate_trial(): a * (x^z)
        escala_gamma = 1.0 / (lambda * k);                    // rgamma(rate = lambda * k)
    }

    // Parametros de la lognormal (una vez por llamada)
    const bool manip_fija = (s == 0);
    double ln_mu = 0.0, ln_sd = 0.0;
    if (!manip_fija) {
        if (modo == 0) {
            if (!(h > 0)) stop("h debe ser > 0 si s > 0");
            ln_sd = std::sqrt(std::log(1 + (s / h) * (s / h)));   // rlnorm0_cpp
            ln_mu = std::log(h) - 0.5 * ln_sd * ln_sd;
        } else {
            ln_sd = s;                                             // simulate_trial()
            ln_mu = std::log(h) - 0.5 * s * s;
        }
    }

    for (int sim = 0; sim < nsim; ++sim) {
        std::fill(tocado.begin(), tocado.end(), 0);
        int distintos = 0;
        double t = 0.0;

        while (t < T) {
            t += R::rgamma(k, escala_gamma);
            if (modo == 0) { if (t > T) break; }
            else           { if (t >= T) break; }

            int pid;
            if (modo == 0) pid = (int) std::floor(R::runif(0, H));
            else           pid = (int) R_unif_index((double) H);

            if (!tocado[pid]) { tocado[pid] = 1; ++distintos; }
            if (corte_saturacion && distintos == H) break;

            double th = manip_fija ? h : std::exp(R::rnorm(ln_mu, ln_sd));
            t += th;
        }

        out[sim] = distintos;
        if (sim % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return out;
}

//' Igual que motor_okuyama_cpp() pero ademas devuelve la cantidad
//' total de vueltas del while (para medir la trampa computacional).
//' Solo para diagnostico (R/00b_motor_rcpp_pruebas.R).
// [[Rcpp::export]]
double contar_vueltas_cpp(int H, double a, double h, double z,
                          double k, double s, double T, int nsim,
                          bool corte_saturacion, double max_vueltas = 1e9) {
    RNGScope scope;
    double rate = k * a * std::pow((double) H, z);
    double escala_gamma = 1.0 / rate;
    const bool manip_fija = (s == 0);
    double ln_sd = manip_fija ? 0 : std::sqrt(std::log(1 + (s / h) * (s / h)));
    double ln_mu = manip_fija ? 0 : std::log(h) - 0.5 * ln_sd * ln_sd;
    std::vector<unsigned char> tocado(H);
    double vueltas = 0;
    for (int sim = 0; sim < nsim; ++sim) {
        std::fill(tocado.begin(), tocado.end(), 0);
        int distintos = 0; double t = 0.0;
        while (t < T) {
            t += R::rgamma(k, escala_gamma);
            if (t > T) break;
            int pid = (int) std::floor(R::runif(0, H));
            if (!tocado[pid]) { tocado[pid] = 1; ++distintos; }
            vueltas += 1;
            if (vueltas >= max_vueltas) return R_PosInf;   // "se colgaria"
            if (corte_saturacion && distintos == H) break;
            t += manip_fija ? h : std::exp(R::rnorm(ln_mu, ln_sd));
        }
        if (sim % 1000 == 0) Rcpp::checkUserInterrupt();
    }
    return vueltas;
}
