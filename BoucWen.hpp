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
    template <typename precision_t>
    class BoucWen_BaseClass {
    public:
      virtual ~BoucWen_BaseClass();
      virtual void List() const = 0;
      virtual precision_t Solver(const precision_t Value) = 0;
    };

    template <typename precision_t>
    BoucWen_BaseClass<precision_t>::~BoucWen_BaseClass(){}

    //=======================================================//
    // Bouc-Wen without degradation                          //
    //=======================================================//
    template <typename precision_t>
    class BoucWen : public BoucWen_BaseClass<precision_t> {
    private:
      precision_t z_old, z_new_p, z_new, z_eval;
      precision_t e_new, e_new_, e_old;
      precision_t Phi, Phi_, Psi;
      precision_t sign;
      precision_t fz_new, fz_new_;
	       
      precision_t tolerance, startPoint;
      int32_t maxIter, count;

      precision_t DispT, deltaDisp;
      precision_t force;
    protected:
      precision_t alpha;    /*!< \brief Post-yield stiffness ratio \f$\alpha = k_y/k_e\f$ with \f$k_y\$ the post
			     * yield stiffness and \f$ke\f$ de pre-yield stiffness. */
      precision_t Fy;       /*!< \brief Yield force */
      precision_t ko;       /* ko is the elasttic stiffness \f$k_0 = F_y/u_y\f$ where \f$F_y\f$ is the post
			     * yield stiffness and \f$u_y\f$ the yield displacement. */
      precision_t beta;     /*!< \brief Bouc-Wen model coefficient. */
      precision_t gamma;    /*!< \brief Bouc-Wen model coefficient. */
      precision_t n;        /*!< \brief Hardening - Softening parameter. Controls the transition from linear to
			     * non-linear range (as n increases the transition becomes sharper; n is usually
			     * grater or equal to 1). */
    public:	       
      BoucWen();
      BoucWen(const precision_t *const parameters);
      ~BoucWen();
      void List() const override;
      precision_t Solver(const precision_t Value) override;
      precision_t signum(const precision_t Value);
    };
	  
    template <typename precision_t>
    BoucWen<precision_t>::BoucWen(const precision_t *const parameters) {
      alpha = parameters[0];
      Fy = parameters[1];
      ko = parameters[2];
      beta = parameters[3];
      gamma = parameters[4];
      n = parameters[5];

      maxIter = 100;
      tolerance = static_cast<precision_t>(1E-12);
      startPoint = static_cast<precision_t>(0.01);
      DispT = static_cast<precision_t>(0.0);
      e_old = static_cast<precision_t>(0.0);
      z_old = static_cast<precision_t>(0.0);
    }
	  
    template <typename precision_t>
    BoucWen<precision_t>::~BoucWen() {}

    template <typename precision_t>
    void BoucWen<precision_t>::List() const {
      std::cout << "Alpha: " << alpha << " ";
      std::cout << "Yield force: " << Fy << " ";
      std::cout << "Elastic stiffness: " << ko << " ";
      std::cout << "Beta: " << beta << " ";
      std::cout << "Gamma: " << gamma << " ";
      std::cout << "n: " << n << " " << std::endl;
    }
  
    template <typename precision_t>
    precision_t BoucWen<precision_t>::signum(const precision_t Value) {
      if (Value > static_cast<precision_t>(0.0)){
	return static_cast<precision_t>(1.0);
      } else if (Value < static_cast<precision_t>(0.0)){
	return (precision_t) -1.0;
      } else return static_cast<precision_t>(0.0);
    }
	  
    template <typename precision_t>
    precision_t BoucWen<precision_t>::Solver(const precision_t DispTdT) {
      deltaDisp = DispTdT - DispT;
	       
      /* Perform Newton-Rhapson */
      count = 0;
      z_new = static_cast<precision_t>(1.0);
      z_new_p = startPoint; z_eval = startPoint;
	       
      while ((std::fabs(z_new_p - z_new) > tolerance) && (count < maxIter)) {
	e_new = e_old + (static_cast<precision_t>(1.0) - alpha)*deltaDisp*ko/Fy*z_eval;

	sign = signum(deltaDisp*z_eval);
	Psi = gamma + beta*sign;
	Phi = static_cast<precision_t>(1.0) - std::pow(std::fabs(z_eval),n)*Psi;

	fz_new = z_eval - z_old - Phi*deltaDisp*ko/Fy;


	/* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
	e_new_ = (static_cast<precision_t>(1.0) - alpha)*ko/Fy*deltaDisp;
	sign = signum(z_eval);
	Phi_ = -n*std::pow(std::fabs(z_eval), n - static_cast<precision_t>(1.0))*sign*Psi;
	fz_new_ = static_cast<precision_t>(1.0) - Phi_*deltaDisp*ko/Fy;

		       
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
	force = alpha*ko*DispTdT + (static_cast<precision_t>(1.0) - alpha)*Fy*z_eval;
		    
	// Compute material degradation parameters.
	e_new = e_old + (static_cast<precision_t>(1.0) - alpha)*ko/Fy*deltaDisp*z_eval;
      }

      DispT = DispTdT;
      e_old = e_new;
      z_old = z_eval;
	       
      return force;
    }

    //=======================================================//
    // Bouc-Wen with degradation                             //
    //=======================================================//
    template <typename precision_t>
    class BoucWen_Deg : public BoucWen<precision_t> {
    private:
      /* Declaration of routine variables */
      precision_t z_old, z_new_p, z_new, z_eval;
      precision_t e_new, e_new_, e_old;
      precision_t A_new, A_new_, nu_new, nu_new_, eta_new, eta_new_;
      precision_t Phi, Phi_, Psi;
      precision_t sign;
      precision_t fz_new, fz_new_;
     
      precision_t tolerance, startPoint;
      int32_t maxIter, count;

      precision_t DispT, deltaDisp;
      precision_t force;
    protected:
      precision_t A0;       /*!< \brief Hysteresis amplitude. */
      precision_t deltaA;   /*!< \brief Control parameter of the hysteresis amplitude with respect to the
			     * energy. */
      precision_t nu0;      /*!< \brief Strength degradation. */
      precision_t deltaNu;  /*!< \brief Strength degradation parameter. With \f$\delta_\nu = 0\f$ no strength
			     * degradation is included in the model. */
      precision_t eta0;     /*!< \brief Stiffness degradation. */
      precision_t deltaEta; /*!< \brief Stiffness degradation parameter. With \f$\delta_\eta = 0\f$ no stiffness
			     * degradation is included in the model.*/
    public:
      BoucWen_Deg();
      ~BoucWen_Deg();
      BoucWen_Deg(precision_t *const parameters);
      precision_t Solver(const precision_t Value) override;
      void List() const override;
    };

    template <typename precision_t>
    BoucWen_Deg<precision_t>::BoucWen_Deg() : BoucWen<precision_t>() {}

    template <typename precision_t>
    BoucWen_Deg<precision_t>::~BoucWen_Deg() {}
	  
    template <typename precision_t>
    BoucWen_Deg<precision_t>::BoucWen_Deg(precision_t *const parameters) : BoucWen<precision_t>(parameters) {
      A0 = parameters[6];
      deltaA = parameters[7];
      nu0 = parameters[8];
      deltaNu = parameters[9];
      eta0 = parameters[10];
      deltaEta = parameters[11];

      maxIter = 100;
      tolerance = static_cast<precision_t>(1E-12);
      startPoint = static_cast<precision_t>(0.01);     
      DispT = static_cast<precision_t>(0.0);
      e_old = static_cast<precision_t>(0.0);
      z_old = static_cast<precision_t>(0.0);       
    }

    template <typename precision_t>
    precision_t BoucWen_Deg<precision_t>::Solver(const precision_t DispTdT) {
      deltaDisp = DispTdT - DispT;
	       
      /* Perform Newton-Rhapson */
      count = 0;
      z_new = 1.0; z_new_p = startPoint; z_eval = startPoint;
	  
      while ( fabs(z_new_p - z_new) > tolerance && count < maxIter ){
	e_new = e_old + (static_cast<precision_t>(1.0) - this->alpha)*deltaDisp*this->ko/this->Fy*z_eval;

	A_new = A0 - deltaA*e_new;
	nu_new = nu0 + deltaNu*e_new;
	eta_new = eta0 + deltaEta*e_new;
 
	sign = this->signum(deltaDisp*z_eval);
	Psi = this->gamma + this->beta*sign;
	Phi = A_new - pow(fabs(z_eval),this->n)*Psi*nu_new;

	fz_new = z_eval - z_old - Phi/eta_new*deltaDisp*this->ko/this->Fy;

	/* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
	e_new_ = (static_cast<precision_t>(1.0) - this->alpha)*this->ko/this->Fy*deltaDisp;
	A_new_ = -deltaA*e_new_;
	nu_new_ = deltaNu*e_new_;
	eta_new_ = deltaEta*e_new_;
	sign = this->signum(z_eval);
	Phi_ = A_new_ - this->n*pow(fabs(z_eval), this->n - static_cast<precision_t>(1.0))*sign*Psi*nu_new - pow(fabs(z_eval), this->n)*Psi*nu_new_;
	fz_new_ = static_cast<precision_t>(1.0) - (Phi_*eta_new - Phi*eta_new_)/pow(eta_new, static_cast<precision_t>(2.0))*deltaDisp*this->ko/this->Fy;

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
	force = this->alpha*this->ko*DispTdT + (static_cast<precision_t>(1.0) - this->alpha)*this->Fy*z_eval;

	// Compute material degradation parameters.
	e_new = e_old + (static_cast<precision_t>(1.0) - this->alpha)*this->ko/this->Fy*deltaDisp*z_eval;
	A_new = A0 - deltaA*e_new;
	nu_new = nu0 + deltaNu*e_new;
	eta_new = eta0 + deltaEta*e_new;
      }

      DispT = DispTdT;
      e_old = e_new;
      z_old = z_eval;

      return force;
    }

    template <typename precision_t>
    void BoucWen_Deg<precision_t>::List() const {
      BoucWen<precision_t>::List();

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
    template <typename precision_t>
    class BoucWen_BaberNoori : public BoucWen_Deg<precision_t> {
    private:
      precision_t z_old, z_new_p, z_new, z_eval;
      precision_t e_new, e_new_, e_old;
      precision_t A_new, A_new_, nu_new, nu_new_, eta_new, eta_new_;
      precision_t Phi, Phi_, Theta, hze, hze_exp, vs_1, vs_2;
      precision_t hze_exp_, vs_1_, hze_1, hze_2, hze_, zu;
      precision_t sign;
      precision_t fz_new, fz_new_;
     
      precision_t tolerance, startPoint;
      int32_t maxIter, count;

      precision_t DispT, deltaDisp;
      precision_t force;
    protected:
      precision_t vs0;      /*!< \brief Pinching severity. With \f$\zeta_s = 0\f$ there is no pinching effect
			     * included in the model. */
      precision_t p;        /*!< \brief Initial pinching parameter. With \f$p = 0\f$ there is no pinching effect
			     * included in the model. */
      precision_t q;        /*!< \brief Pinching parameter. */
      precision_t psi0;     /*!< \brief Pinching parameter. */
      precision_t deltaPsi; /*!< \brief Controls the change of pinching in the model */
      precision_t lambda;   /*!< \brief Pinching parameter. */
    public:
      BoucWen_BaberNoori();
      ~BoucWen_BaberNoori();
      BoucWen_BaberNoori(precision_t *const parameters);
      precision_t Solver(const precision_t Value) final;
      void List() const final;
    };

    template <typename precision_t>
    BoucWen_BaberNoori<precision_t>::BoucWen_BaberNoori() : BoucWen_Deg<precision_t>() {}

    template <typename precision_t>
    BoucWen_BaberNoori<precision_t>::~BoucWen_BaberNoori() {}
	  
    template <typename precision_t>
    BoucWen_BaberNoori<precision_t>::BoucWen_BaberNoori(precision_t *const parameters) : BoucWen_Deg<precision_t>(parameters) {
      vs0 = parameters[12];
      p = parameters[13];
      q = parameters[14];
      psi0 = parameters[15];
      deltaPsi = parameters[16];
      lambda = parameters[17];

      maxIter = 100;
      tolerance = static_cast<precision_t>(1E-12);
      startPoint = static_cast<precision_t>(0.01);
      DispT = static_cast<precision_t>(0.0);
      e_old = static_cast<precision_t>(0.0);
      z_old = static_cast<precision_t>(0.0); 
    }

    template <typename precision_t>
    precision_t BoucWen_BaberNoori<precision_t>::Solver(const precision_t DispTdT) {
      deltaDisp = DispTdT - DispT;

      /* Perform Newton-Rhapson */
      count = 0;
      z_new = 1.0; z_new_p = startPoint; z_eval = startPoint;
	  
      while ( fabs(z_new_p - z_new) > tolerance && count < maxIter ){
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
	hze_exp = exp(-std::pow(z_eval*this->signum(deltaDisp) -q*zu, 2.0)/std::pow(vs_2, 2.0));
	hze = 1.0 - vs_1*hze_exp;
	       
	fz_new = z_eval - z_old - hze*Phi/eta_new*deltaDisp;

	/* Evaluate function derivatives with respect to z_eval for the Newton-Rhapson scheme */
	e_new_ = (1.0 - this->alpha)*this->ko/this->Fy*deltaDisp;
	A_new_ = -this->deltaA*e_new_;
	nu_new_ = this->deltaNu*e_new_;
	eta_new_ = this->deltaEta*e_new_;
	sign = this->signum(z_eval);
	Phi_ = A_new_ - this->n*std::pow(fabs(z_eval), this->n - 1.0)*sign*Theta*nu_new - std::pow(fabs(z_eval), this->n)*Theta*nu_new_;
	       
	hze_ = deltaDisp*this->ko/this->Fy*p*vs0*exp(-std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))), 1/this->n),2)/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))*(this->alpha - 1) - vs0*exp(-std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2)/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1)*((2*(this->signum(deltaDisp) - (this->deltaNu*deltaDisp*this->ko/this->Fy*q*(this->alpha - 1)*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n - 1))/(this->n*(this->beta + this->gamma)*std::pow(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)))*(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n)))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)) + (2*deltaPsi*deltaDisp*this->ko/this->Fy*(this->alpha - 1)*std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),2)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),3)) + (2*deltaDisp*this->ko/this->Fy*p*vs0*exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))*(this->alpha - 1)*std::pow(z_eval*this->signum(deltaDisp) - q*std::pow(1/((this->beta + this->gamma)*(this->nu0 + this->deltaNu*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)))),1/this->n),2))/(std::pow(lambda - vs0*(exp(-p*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1))) - 1),3)*std::pow(psi0 + deltaPsi*(e_old - deltaDisp*this->ko/this->Fy*z_eval*(this->alpha - 1)),2)));
	       
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
	  
    template <typename precision_t>
    void BoucWen_BaberNoori<precision_t>::List() const {
      BoucWen_Deg<precision_t>::List();

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
