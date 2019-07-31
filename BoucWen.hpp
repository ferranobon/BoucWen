#ifndef _BOUCWEN_HPP_
#define _BOUCWEN_HPP_

#include <iostream>
#include <cmath>
#include <cstdint>

namespace BoucWen {
    namespace NewtonRhapson {

        //=======================================================//
        // Base class definition                                 //
        //=======================================================//
        template <typename T>
        class BoucWen_BaseClass {
        public:
            virtual ~BoucWen_BaseClass();
            virtual void List() const = 0;
            virtual T Solver(const T Value) = 0;
        };

        template <typename T>
        BoucWen_BaseClass<T>::~BoucWen_BaseClass(){}

        //=======================================================//
        // Bouc-Wen without degradation                          //
        //=======================================================//
        template <typename T>
        class BoucWen : public BoucWen_BaseClass<T> {
        private:
            T z_old, z_new_p, z_new, z_eval;
            T e_new, e_new_, e_old;
            T Phi, Phi_, Psi;
            T sign;
            T fz_new, fz_new_;

            T tolerance, startPoint;
            int32_t maxIter, count;

            T DispT, deltaDisp;
            T force;
        protected:
            T alpha;    /*!< \brief Post-yield stiffness ratio \f$\alpha = k_y/k_e\f$ with \f$k_y\$ the post
             * yield stiffness and \f$ke\f$ de pre-yield stiffness. */
            T Fy;       /*!< \brief Yield force */
            T ko;       /* ko is the elasttic stiffness \f$k_0 = F_y/u_y\f$ where \f$F_y\f$ is the post
             * yield stiffness and \f$u_y\f$ the yield displacement. */
            T beta;     /*!< \brief Bouc-Wen model coefficient. */
            T gamma;    /*!< \brief Bouc-Wen model coefficient. */
            T n;        /*!< \brief Hardening - Softening parameter. Controls the transition from linear to
             * non-linear range (as n increases the transition becomes sharper; n is usually
             * grater or equal to 1). */
        public:
            BoucWen();
            BoucWen(const T *const parameters);
            ~BoucWen();

            static constexpr uint8_t NUM_PARAMS = 6;

            void List() const override;
            T Solver(const T Value) override;
            T signum(const T Value);
        };

        template <typename T>
        BoucWen<T>::BoucWen(const T *const parameters) {
            alpha = parameters[0];
            Fy = parameters[1];
            ko = parameters[2];
            beta = parameters[3];
            gamma = parameters[4];
            n = parameters[5];

            maxIter = 100;
            tolerance = static_cast<T>(1E-12);
            startPoint = static_cast<T>(0.01);
            DispT = static_cast<T>(0.0);
            e_old = static_cast<T>(0.0);
            z_old = static_cast<T>(0.0);
        }

        template <typename T>
        BoucWen<T>::~BoucWen() {}

        template <typename T>
        void BoucWen<T>::List() const {
            std::cout << "Alpha: " << alpha << " ";
            std::cout << "Yield force: " << Fy << " ";
            std::cout << "Elastic stiffness: " << ko << " ";
            std::cout << "Beta: " << beta << " ";
            std::cout << "Gamma: " << gamma << " ";
            std::cout << "n: " << n << " " << std::endl;
        }

        template <typename T>
        T BoucWen<T>::signum(const T Value) {
            if (Value > static_cast<T>(0.0)){
                return static_cast<T>(1.0);
            } else if (Value < static_cast<T>(0.0)){
                return (T) -1.0;
            } else return static_cast<T>(0.0);
        }

        template <typename T>
        T BoucWen<T>::Solver(const T DispTdT) {
            deltaDisp = DispTdT - DispT;

            /* Perform Newton-Rhapson */
            count = 0;
            z_new = static_cast<T>(1.0);
            z_new_p = startPoint; z_eval = startPoint;

            while ((std::fabs(z_new_p - z_new) > tolerance) && (count < maxIter)) {
                e_new = e_old + (static_cast<T>(1.0) - alpha)*deltaDisp*ko/Fy*z_eval;

                sign = signum(deltaDisp*z_eval);
                Psi = gamma + beta*sign;
                Phi = static_cast<T>(1.0) - std::pow(std::fabs(z_eval),n)*Psi;

                fz_new = z_eval - z_old - Phi*deltaDisp*ko/Fy;


                /* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
                e_new_ = (static_cast<T>(1.0) - alpha)*ko/Fy*deltaDisp;
                sign = signum(z_eval);
                Phi_ = -n*std::pow(std::fabs(z_eval), n - static_cast<T>(1.0))*sign*Psi;
                fz_new_ = static_cast<T>(1.0) - Phi_*deltaDisp*ko/Fy;


                /* Perform a new step */
                z_new = z_eval - fz_new/fz_new_;

                /* Update the root */
                z_new_p = z_eval;
                z_eval = z_new;

                count = count + 1;

                /* Warning if there is no convergence */
                if (count == maxIter) {
                    std::cerr << "WARNING: BoucWen() -- could not find the root z_{i+1}, after " << maxIter << " iterations and norm: " << std::fabs(z_new_p - z_new) << std::endl;
                }

                // Compute restoring force.
                force = alpha*ko*DispTdT + (static_cast<T>(1.0) - alpha)*Fy*z_eval;

                // Compute material degradation parameters.
                e_new = e_old + (static_cast<T>(1.0) - alpha)*ko/Fy*deltaDisp*z_eval;
            }

            DispT = DispTdT;
            e_old = e_new;
            z_old = z_eval;

            return force;
        }

        //=======================================================//
        // Bouc-Wen with degradation                             //
        //=======================================================//
        template <typename T>
        class BoucWen_Deg : public BoucWen<T> {
        private:
            /* Declaration of routine variables */
            T z_old, z_new_p, z_new, z_eval;
            T e_new, e_new_, e_old;
            T A_new, A_new_, nu_new, nu_new_, eta_new, eta_new_;
            T Phi, Phi_, Psi;
            T sign;
            T fz_new, fz_new_;

            T tolerance, startPoint;
            int32_t maxIter, count;

            T DispT, deltaDisp;
            T force;
        protected:
            T A0;       /*!< \brief Hysteresis amplitude. */
            T deltaA;   /*!< \brief Control parameter of the hysteresis amplitude with respect to the
                         * energy. */
            T nu0;      /*!< \brief Strength degradation. */
            T deltaNu;  /*!< \brief Strength degradation parameter. With \f$\delta_\nu = 0\f$ no strength
                         * degradation is included in the model. */
            T eta0;     /*!< \brief Stiffness degradation. */
            T deltaEta; /*!< \brief Stiffness degradation parameter. With \f$\delta_\eta = 0\f$ no stiffness
                         * degradation is included in the model.*/
        public:
            BoucWen_Deg();
            ~BoucWen_Deg();
            BoucWen_Deg(T *const parameters);

            static constexpr uint8_t NUM_PARAMS = 12;

            T Solver(const T Value) override;
            void List() const override;
        };

        template <typename T>
        BoucWen_Deg<T>::BoucWen_Deg() : BoucWen<T>() {}

        template <typename T>
        BoucWen_Deg<T>::~BoucWen_Deg() {}

        template <typename T>
        BoucWen_Deg<T>::BoucWen_Deg(T *const parameters) : BoucWen<T>(parameters) {
            A0 = parameters[6];
            deltaA = parameters[7];
            nu0 = parameters[8];
            deltaNu = parameters[9];
            eta0 = parameters[10];
            deltaEta = parameters[11];

            maxIter = 100;
            tolerance = static_cast<T>(1E-12);
            startPoint = static_cast<T>(0.01);
            DispT = static_cast<T>(0.0);
            e_old = static_cast<T>(0.0);
            z_old = static_cast<T>(0.0);
        }

        template <typename T>
        T BoucWen_Deg<T>::Solver(const T DispTdT) {
            deltaDisp = DispTdT - DispT;

            /* Perform Newton-Rhapson */
            count = 0;
            z_new = 1.0; z_new_p = startPoint; z_eval = startPoint;

            while ( fabs(z_new_p - z_new) > tolerance && count < maxIter ){
                e_new = e_old + (static_cast<T>(1.0) - this->alpha)*deltaDisp*this->ko/this->Fy*z_eval;

                A_new = A0 - deltaA*e_new;
                nu_new = nu0 + deltaNu*e_new;
                eta_new = eta0 + deltaEta*e_new;

                sign = this->signum(deltaDisp*z_eval);
                Psi = this->gamma + this->beta*sign;
                Phi = A_new - pow(fabs(z_eval),this->n)*Psi*nu_new;

                fz_new = z_eval - z_old - Phi/eta_new*deltaDisp*this->ko/this->Fy;

                /* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
                e_new_ = (static_cast<T>(1.0) - this->alpha)*this->ko/this->Fy*deltaDisp;
                A_new_ = -deltaA*e_new_;
                nu_new_ = deltaNu*e_new_;
                eta_new_ = deltaEta*e_new_;
                sign = this->signum(z_eval);
                Phi_ = A_new_ - this->n*pow(fabs(z_eval), this->n - static_cast<T>(1.0))*sign*Psi*nu_new - pow(fabs(z_eval), this->n)*Psi*nu_new_;
                fz_new_ = static_cast<T>(1.0) - (Phi_*eta_new - Phi*eta_new_)/pow(eta_new, static_cast<T>(2.0))*deltaDisp*this->ko/this->Fy;

                /* Perform a new step */
                z_new = z_eval - fz_new/fz_new_;

                /* Update the root */
                z_new_p = z_eval;
                z_eval = z_new;

                count = count + 1;

                /* Warning if there is no convergence */
                if (count == maxIter) {
                    std::cerr << "WARNING: BoucWen() -- could not find the root z_{i+1}, after " << maxIter << " iterations and norm: " << std::fabs(z_new_p - z_new) << std::endl;
                }

                // Compute restoring force.
                force = this->alpha*this->ko*DispTdT + (static_cast<T>(1.0) - this->alpha)*this->Fy*z_eval;

                // Compute material degradation parameters.
                e_new = e_old + (static_cast<T>(1.0) - this->alpha)*this->ko/this->Fy*deltaDisp*z_eval;
                A_new = A0 - deltaA*e_new;
                nu_new = nu0 + deltaNu*e_new;
                eta_new = eta0 + deltaEta*e_new;
            }

            DispT = DispTdT;
            e_old = e_new;
            z_old = z_eval;

            return force;
        }

        template <typename T>
        void BoucWen_Deg<T>::List() const {
            BoucWen<T>::List();

            std::cout << "A0: " << A0 << " ";
            std::cout << "deltaA: " << deltaA << " ";
            std::cout << "nu0: " << nu0 << " ";
            std::cout << "deltaNu: " << deltaNu << " ";
            std::cout << "eta0: " << eta0 << " ";
            std::cout << "deltaEta: " << deltaEta << " " << std::endl;
        }

        //=======================================================//
        // Bouc-Wen with pitching and degradation (Baber-Noori)  //
        //=======================================================//
        template <typename T>
        class BoucWen_BaberNoori : public BoucWen_Deg<T> {
        private:
            T z_old, z_new_p, z_new, z_eval;
            T e_new, e_new_, e_old;
            T A_new, A_new_, nu_new, nu_new_, eta_new, eta_new_;
            T Phi, Phi_, Theta, hze, hze_exp, vs_1, vs_2;
            T hze_exp_, vs_1_, hze_1, hze_2, hze_, zu;
            T sign;
            T fz_new, fz_new_;

            T tolerance, startPoint;
            int32_t maxIter, count;

            T DispT, deltaDisp;
            T force;
        protected:
            T vs0;      /*!< \brief Pinching severity. With \f$\zeta_s = 0\f$ there is no pinching effect
             * included in the model. */
            T p;        /*!< \brief Initial pinching parameter. With \f$p = 0\f$ there is no pinching effect
             * included in the model. */
            T q;        /*!< \brief Pinching parameter. */
            T psi0;     /*!< \brief Pinching parameter. */
            T deltaPsi; /*!< \brief Controls the change of pinching in the model */
            T lambda;   /*!< \brief Pinching parameter. */
        public:
            BoucWen_BaberNoori();
            ~BoucWen_BaberNoori();
            BoucWen_BaberNoori(T *const parameters);

            static constexpr uint8_t NUM_PARAMS = 18;

            T Solver(const T Value) final;
            void List() const final;
        };

        template <typename T>
        BoucWen_BaberNoori<T>::BoucWen_BaberNoori() : BoucWen_Deg<T>() {}

        template <typename T>
        BoucWen_BaberNoori<T>::~BoucWen_BaberNoori() {}

        template <typename T>
        BoucWen_BaberNoori<T>::BoucWen_BaberNoori(T *const parameters) : BoucWen_Deg<T>(parameters) {
            vs0 = parameters[12];
            p = parameters[13];
            q = parameters[14];
            psi0 = parameters[15];
            deltaPsi = parameters[16];
            lambda = parameters[17];

            maxIter = 100;
            tolerance = static_cast<T>(1E-12);
            startPoint = static_cast<T>(0.01);
            DispT = static_cast<T>(0.0);
            e_old = static_cast<T>(0.0);
            z_old = static_cast<T>(0.0);
        }

        template <typename T>
        T BoucWen_BaberNoori<T>::Solver(const T DispTdT) {
            deltaDisp = DispTdT - DispT;

            /* Perform Newton-Rhapson */
            count = 0;
            z_new = 1.0; z_new_p = startPoint; z_eval = startPoint;

            while(std::fabs(z_new_p - z_new) > tolerance && count < maxIter){
                e_new = e_old + (1.0 - this->alpha)*this->ko/this->Fy*deltaDisp*z_eval;
                A_new = this->A0 - this->deltaA*e_new;
                nu_new = this->nu0 + this->deltaNu*e_new;
                eta_new = this->eta0 + this->deltaEta*e_new;

                sign = this->signum(deltaDisp*z_eval);
                Theta = this->gamma + this->beta*sign;
                Phi = A_new - std::pow(fabs(z_eval),this->n)*Theta*nu_new;

                vs_1 = vs0*(1.0 - exp(-p*e_new));
                vs_2 = (psi0 + deltaPsi*e_new)*(lambda + vs_1);
                zu = std::pow(1.0/((this->nu0 + this->deltaNu*e_new)*(this->beta + this->gamma)), 1.0/this->n);
                hze_exp = std::exp(-std::pow(z_eval*this->signum(deltaDisp) -q*zu, 2.0)/std::pow(vs_2, 2.0));
                hze = 1.0 - vs_1*hze_exp;

                fz_new = z_eval - z_old - hze*Phi/eta_new*deltaDisp;

                /* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
                e_new_ = (1.0 - this->alpha)*this->ko/this->Fy*deltaDisp;
                A_new_ = -this->deltaA*e_new_;
                nu_new_ = this->deltaNu*e_new_;
                eta_new_ = this->deltaEta*e_new_;
                sign = this->signum(z_eval);
                Phi_ = A_new_ - this->n*std::pow(std::fabs(z_eval), this->n - 1.0)*sign*Theta*nu_new - std::pow(std::fabs(z_eval), this->n)*Theta*nu_new_;

                hze_ = deltaDisp*this->ko/this->Fy*p*vs0*std::exp(-std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))), 1/this->n),2)/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))*(this->alpha - 1) - vs0*exp(-std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2)/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1)*((2*(this->signum(deltaDisp) - (this->deltaNu*deltaDisp*this->ko/this->Fy*q*(this->alpha - 1)*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n - 1))/(this->n*(this->beta + this->gamma)*std::pow(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n)))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)) + (2*deltaPsi*deltaDisp*this->ko/this->Fy*(this->alpha - 1)*std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),3)) + (2*deltaDisp*this->ko/this->Fy*p*vs0*exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))*(this->alpha - 1)*std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),3)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)));

                fz_new_ = 1.0 - (hze_*Phi*eta_new + hze*Phi_*eta_new - hze*Phi*eta_new_)/std::pow(eta_new, 2.0)*deltaDisp;

                /* Perform a new step */
                z_new = z_eval - fz_new/fz_new_;

                /* Update the root */
                z_new_p = z_eval;
                z_eval = z_new;

                count = count + 1;

                /* Warning if there is no convergence */
                if (count == maxIter) {
                    std::cerr << "WARNING: BoucWen() -- could not find the root z_{i+1}, after " << maxIter << " iterations and norm: " << std::fabs(z_new_p - z_new) << std::endl;
                }

                // Compute restoring force.
                force = this->alpha*this->ko*DispTdT + (1.0 - this->alpha)*this->Fy*z_eval;

                // Compute material degradation parameters.
                e_new = e_old + (1.0 - this->alpha)*this->ko/this->Fy*deltaDisp*z_eval;
                A_new = this->A0 - this->deltaA*e_new;
                nu_new = this->nu0 + this->deltaNu*e_new;
                eta_new = this->eta0 + this->deltaEta*e_new;
            }

            DispT = DispTdT;
            e_old = e_new;
            z_old = z_eval;

            return force;
        }

        template <typename T>
        void BoucWen_BaberNoori<T>::List() const {
            BoucWen_Deg<T>::List();

            std::cout << "vs0: " << vs0 << " ";
            std::cout << "p: " << p << " ";
            std::cout << "q: " << q << " ";
            std::cout << "psi00: " << psi0 << " ";
            std::cout << "deltaPsi: " << deltaPsi << " ";
            std::cout << "lambda: " << lambda << " " << std::endl;
        }
    };
};

#endif //_BOUCWEN_HPP_
