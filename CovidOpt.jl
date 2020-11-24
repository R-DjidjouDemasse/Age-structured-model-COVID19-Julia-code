##### This code was implemented by Quentin Richard, Samuel Alizon, Marc Choisy, Mircea T. Sofonea and Ramses Djidjou-Demasse #####
##### The model (and some results) are described in Q. Richard et al (2020) medRxiv: https://www.medrxiv.org/content/10.1101/2020.06.23.20138099v2 #####



############### To execute first ###############

# cd("C:/.../Julia") # to access the work directory

#### Packages needed ####

# using Plots; pyplot()	# to draw curves
# using LinearAlgebra 	# to use linear algebra
# using Plots.PlotMeasures 	# (to put margins in cm, ex :left_margin=10cm)))

#### To load the file ####

# include("CovidOpt.jl")

############### Functions used in the following ###############

function Charac(y,a,b)	# indicator function on [a,b]
	if b>y>=a
		return 1
	else
		return 0
	end
end

function Expo(y,gamma)	# decreasing exponential function
	exp.(-y*gamma);
end

function Gauss(y,mu,sigma)	# Gaussiafn distribution with mean mu and variance sigma
	(1/(sigma*sqrt(2*pi)))*exp.(-0.5*(((y-mu*ones(size(y)))./sigma).^2));
end

function Nroot(n,A)	# function N root (with n>1), i.e. x^N=y <=> x=Nroot(y)
	x=1; 			# initial condition
	eps=10^(-6); 	# accuracy
	err=abs(x^n-A);	# error

	while err>eps
		x=(1/n)*((n-1)*x+(A/(x^(n-1))));
		err=abs(x^n-A);
	end
	return x
end

function Gamma(x) # Gamma distribution
	t=collect(0:0.001:10000);
	Nt=length(t);
	func=(t.^(x-1)).*exp.(-t);
	Integ=zeros(Nt-1);
	for i=1:(Nt-1)
		Integ[i]=0.0005*(func[i]+func[i+1]);
	end
	Sum=sum(Integ);
	return Sum
end

############### The PDE model ###############

	##### Algorithm and figures mortality rates, gamma and mu (Figure 2.b) #####

function Covid(T,c,Hsat,p,isympt,xip,country)

	## Inputs: T: maximal time (in days), c: control measure (put 0 for no control), p: proportion of paucisymptomatic, isympt: mean time to develop symptom
	## Country="France", "Burkina" or "Vietnam"
	## xip : reduction rate of transmission of COVID-19 in paucisymptomatic population
	## Hsat: ICU hospital beds (5000 for France, 11 for Burkina, 5932 for Vietnam)

		# Time;

	dt=0.2;				# discretisation time (in days)
	nbr=Int(T/dt)+1;	# number of iterations
	t=collect(0:dt:T);	# vector of time

		# Age

	amin=1;		# discretisation in age

	if country=="France"
		amax=20;	# (20 groups of 5 years age, first one is 0-5 years, 2nd: 5-10, etc. - last one is 95 and more)
	elseif country=="Burkina"
		amax=18;	# (18 groups of 5 year age, - last one is 85 and more-)
	elseif country=="Vietnam"
		amax=20;	# (20 groups of 5 years age, first one is 0-5 years, 2nd: 5-10, etc. - last one is 95 and more)
	else
		println("You did not write the country correctly, type either France, Burkina, or Vietnam")
		return
	end

	Na=Int(amax); 			# number of groups (20)
	a=collect(amin:1:amax); # vector of age

	if c==0
		c=zeros(nbr,Na);	# if we put no control (0) then we create a matrix of 0;
	end

		# Age of infection

	imin=0; imax=34; imaxS=25.2;  imaxM=22.2; 	# min and max of age of infection (S for severe, M for mild)
	di=0.2;             						# discretisation age of infection
	Ni=Int(imax*(1/di))+1;	 				# number of iterations
	NiS=Int(imaxS*(1/di))+1;	 			# number of age iterations until end of symptomes in severe infected
	NiM=Int(imaxM*(1/di))+1;	 			# number of age iterations until end of symptomes in mild infected (and paucisymptomatic)
	i=collect(imin:di:imax);				# vector of age of infection
	Tsympt=Int(isympt*(1/di))+1;			# number that corresponds to the appearance of stymptoms
	NiSympt=Int(Int((10*imaxS-10*isympt)*(1/di))/10)+1;   # number of iterations between appearance of symptoms and end of symptom in severe infected

		# Parameters

	munat=zeros(Na); 				# natural mortality rate by day
	muadd=zeros(Na); 				# additional natural mortality (due to saturation of hospital) by day
	gammaIdir=zeros(Na,Ni);			# mortality rate directly due to COVID-19, depending in age and time since infection (of hospitalised patients)
	gammaIindir=zeros(Na,Ni);		# mortality rate indirectly due to COVID-19

	betaS=zeros(1,Ni);				# transmission rate for severe
	betaM=zeros(1,Ni);				# transmission rate for mild
	betaA=zeros(1,Ni);				# transmission rate for paucisymptomatic
	nu=xip;							# reduction rate of transmission of COVID-19 in paucisymptomatic population

	xiS=ones(Ni);					# case isolation of severe
	xiS[Tsympt:end]=exp.(-log(10)*(collect(isympt:di:imax)-isympt*ones(Ni-Tsympt+1)));
	xiM=ones(Ni);					# case isolation of mild
	xiM[Tsympt:end]=exp.(-log(2)*(collect(isympt:di:imax)-isympt*ones(Ni-Tsympt+1)));

	hS=ones(Na,Ni);   				# recovery rate for severe
	hM=ones(Na,Ni);					# recovery rate for mild (for paucisymptomatic)
	hS[:,1:NiS]=zeros(Na,NiS);
	hM[:,1:NiM]=zeros(Na,NiM);

	prophosp=zeros(Na);				# proportion of symptomatic, according to the age, that need hospitalisation

	if country=="France"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273,0.273,0.273]; # (data : Fergusson (2020))
		munat=[1-Nroot(365,1-0.00023),1-Nroot(365,1-0.00007),1-Nroot(365,1-0.00008),1-Nroot(365,1-0.00021),1-Nroot(365,1-0.00040), 1-Nroot(365,1-0.00045),1-Nroot(365,1-0.00057),1-Nroot(365,1-0.00081),1-Nroot(365,1-0.00125),1-Nroot(365,1-0.00204), 1-Nroot(365,1-0.00332), 1-Nroot(365,1-0.0052), 1-Nroot(365,1-0.0078), 1-Nroot(365,1-0.0106), 1-Nroot(365,1-0.0157), 1-Nroot(365,1-0.0217), 1-Nroot(365,1-0.0506),1-Nroot(365,1-0.0686),1-Nroot(365,1-0.151), 1-Nroot(365,1-0.231)];
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245),1-Nroot(NiSympt,1-0.3545),1-Nroot(NiSympt,1-0.3845)]; 		# I had to smooth a little bit the data in the end, because the last data was 80 years old and more
		alpha=1.349*10^(-7);				# to have R0=3.3 (Pasteur paper and reports of ETE team, Montpellier)
	elseif country=="Burkina"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273]; # (data : Fergusson (2020))
		munat=[1-Nroot(365,1-0.0138),1-Nroot(365,1-0.0031),1-Nroot(365,1-0.0017),1-Nroot(365,1-0.0024),1-Nroot(365,1-0.0034), 1-Nroot(365,1-0.0041),1-Nroot(365,1-0.0048),1-Nroot(365,1-0.0059),1-Nroot(365,1-0.0071),1-Nroot(365,1-0.0092), 1-Nroot(365,1-0.0122), 1-Nroot(365,1-0.0174), 1-Nroot(365,1-0.0252), 1-Nroot(365,1-0.0385), 1-Nroot(365,1-0.0578), 1-Nroot(365,1-0.0847), 1-Nroot(365,1-0.1262),1-Nroot(365,1-0.2032)];
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245)];
		alpha=1.34*10^(-7);					# to have R0=3.3
	elseif country=="Vietnam"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273,0.273,0.273]; # (data : Fergusson (2020))
		munat=(1-Nroot(365,1-6.317/1000))*ones(Na);
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245),1-Nroot(NiSympt,1-0.3545),1-Nroot(NiSympt,1-0.3845)]; 		# I had to smooth a little bit the data in the end, because the last data was 80 years old and more
		alpha=0.605*10^(-7);				# to have R0=3.3 (Pasteur paper and reports of ETE team, Montpellier)
	end

	# Let p the proba of dying naturally in 1 year, then 1-p=(1-q)^{365}, with q the proba to die in 1 day
	# gammaI_age is the proba to die, according to the age, from COVID-19 // gammaI_infec is the proba to die, according to the stage of developpement of the infection
	gammaI_infec=zeros(Ni); 			 	  # distribution of deaths due to COVID-19, uniformly in time since infection
	gammaI_infec[Tsympt:NiS]=ones(NiSympt);	  # we can die only if time since infection is between i_sympt and i_max^S
	gammaIdir=gammaI_age*gammaI_infec';

	betak=3;betalambda=5/(log(2)^(1/betak)); # proba of transmission. Weibull(k,lambda) distribution, so that mean=5, variance=5 and standard deviation 1.9 (Ferreti et al-2020)
	betaAA=(alpha)*(betak/betalambda)*((i/betalambda).^(betak-1)).*exp.(-(i/betalambda).^betak);
	betaS=betaAA.*xiS;
    betaM=betaAA.*xiM;
	betaA=betaAA*nu;

	betacontact_age=ones(Na,Na);	# Social matrix contacts, data from Mossong et al (2008) for France, Prem, Cool and Jit (2017) for Burkina and Vietnam


	if country=="France"
		betacontact_age=[3.805	 1.106	 0.412	 0.255	 0.342	 0.753	 1.209	 1.042	 0.496	 0.31	 0.307	 0.272	 0.186	 0.128	 0.088	 0.054	 0.054	 0.054	 0.054	 0.054;
			1.062	 5.033	 1.011	 0.273	 0.164	 0.415	 0.921	 1.138	 0.845	 0.331	 0.207	 0.185	 0.177	 0.112	 0.058	 0.053	 0.053	 0.053	 0.053	 0.053;
			0.238	 1.538	 6.986	 0.846	 0.287	 0.244	 0.403	 0.827	 1.081	 0.57	 0.266	 0.147	 0.102	 0.088	 0.07	 0.07	 0.07	 0.07	 0.07	 0.07;
			0.124	 0.308	 2.301	 7.832	 1.36	 0.651	 0.531	 0.762	 1.025	 1.049	 0.479	 0.164	 0.072	 0.053	 0.035	 0.026	 0.026	 0.026	 0.026	 0.026;
			0.2	 0.168	 0.256	 2.165	 3.934	 1.708	 1.168	 1.133	 0.97	 1.27	 0.813	 0.327	 0.084	 0.046	 0.054	 0.052	 0.052	 0.052	 0.052	 0.052;
			0.533	 0.249	 0.144	 0.731	 1.944	 3.458	 1.777	 1.471	 1.263	 0.997	 0.961	 0.357	 0.113	 0.052	 0.031	 0.024	 0.024	 0.024	 0.024	 0.024;
			0.722	 0.848	 0.53	 0.423	 1.005	 1.72	 2.919	 1.894	 1.457	 1.135	 0.81	 0.405	 0.184	 0.085	 0.047	 0.053	 0.053	 0.053	 0.053	 0.053;
			0.704	 1.072	 0.837	 0.683	 0.765	 1.407	 1.76	 3.211	 2.098	 1.344	 0.92	 0.364	 0.26	 0.152	 0.098	 0.046	 0.046	 0.046	 0.046	 0.046;
			0.31	 0.657	 0.983	 1.109	 0.924	 1.252	 1.591	 1.877	 2.954	 1.682	 1.113	 0.273	 0.184	 0.111	 0.087	 0.048	 0.048	 0.048	 0.048	 0.048;
			0.406	 0.467	 0.612	 1.537	 0.945	 0.994	 1.259	 1.465	 1.573	 2.196	 1.086	 0.341	 0.151	 0.086	 0.083	 0.082	 0.082	 0.082	 0.082	 0.082;
			0.252	 0.595	 0.829	 1.202	 1.019	 1.291	 1.179	 1.215	 1.594	 1.743	 1.918	 0.649	 0.267	 0.119	 0.088	 0.09	 0.09	 0.09	 0.09	 0.09;
			0.542	 0.671	 0.574	 0.749	 0.614	 0.932	 0.918	 0.748	 0.793	 0.681	 0.923	 1.462	 0.513	 0.218	 0.12	 0.095	 0.095	 0.095	 0.095	 0.095;
			0.396	 0.363	 0.257	 0.361	 0.339	 0.508	 0.639	 0.694	 0.523	 0.416	 0.441	 0.675	 1.406	 0.426	 0.263	 0.129	 0.129	 0.129	 0.129	 0.129;
			0.198	 0.313	 0.257	 0.157	 0.212	 0.307	 0.45	 0.454	 0.425	 0.302	 0.324	 0.44	 0.518	 1.003	 0.276	 0.139	 0.139	 0.139	 0.139	 0.139;
			0.105	 0.294	 0.308	 0.355	 0.156	 0.245	 0.253	 0.477	 0.515	 0.434	 0.339	 0.338	 0.66	 0.619	 0.974	 0.326	 0.326	 0.326	 0.326	 0.326;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643];
	elseif country=="Burkina"
		betacontact_age[1:3,:]=[6.415 3.639 1.748 0.919 1.198 1.574 1.469 1.158 0.661 0.398 0.433 0.377 0.251 0.163 0.079 0.036 0.036 0.036;3.477 18.279 4.124 1.077 0.652 1.292 1.3506 1.237 0.911 0.452 0.356 0.345 0.283 0.160 0.068 0.037 0.037 0.037;1.177 5.907 17.935 1.868 0.832 0.717 0.830 0.862 0.851 0.500 0.314 0.169 0.1295 0.113 0.068 0.033 0.033 0.033];
		betacontact_age[4:6,:]=[0.634 1.532 5.808 13.816 2.361 1.144 0.784 0.923 0.852 0.703 0.440 0.228 0.141 0.064 0.031 0.015 0.015 0.015;0.849 0.858 0.911 4.087 4.586 2.153 1.273 1.049 0.824 0.768 0.581 0.346 0.179 0.038 0.029 0.017 0.017 0.017;1.172 0.782 0.525 1.225 2.435 2.846 1.660 1.282 1.038 0.797 0.725 0.407 0.221 0.038 0.013 0.010 0.010 0.010];
		betacontact_age[7:9,:]=[1.125 1.647 1.311 0.664 1.163 1.753 1.842 1.502 1.175 0.902 0.657 0.471 0.248 0.043 0.021 0.013 0.013 0.013;0.900 1.390 1.052 0.789 0.725 1.249 1.356 1.675 1.461 0.962 0.744 0.385 0.200 0.063 0.028 0.009 0.009 0.0091;0.738 1.218 1.200 1.010 0.833 1.094 1.302 1.365 1.516 1.168 0.891 0.366 0.255 0.077 0.029 0.013 0.013 0.013];
		betacontact_age[10:12,:]=[0.464 0.922 0.951 1.188 0.690 0.873 1.009 1.117 1.070 1.035 0.781 0.448 0.186 0.044 0.027 0.022 0.022 0.022;0.565 0.920 1.185 1.005 0.724 1.029 1.013 0.980 1.238 1.241 1.001 0.639 0.234 0.040 0.030 0.022 0.022 0.022;1.290 1.653 1.259 0.956 0.700 1.103 1.183 0.908 0.985 0.923 0.928 0.627 0.346 0.108 0.033 0.038 0.038 0.038];
		betacontact_age[13:15,:]=[0.827 0.941 0.750 0.525 0.443 0.629 0.684 0.775 0.693 0.622 0.575 0.506 0.219 0.098 0.039 0.010 0.010 0.010;0.575 0.935 0.720 0.257 0.232 0.246 0.335 0.357 0.284 0.163 0.189 0.174 0.126 0.083 0.059 0.017 0.017 0.017;0.244 0.792 0.6343 0.393 0.103 0.205 0.141 0.258 0.241 0.238 0.243 0.144 0.154 0.099 0.061 0.061 0.061 0.061];
		betacontact_age[16:18,:]=[0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022;0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022;0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022];
	elseif country=="Vietnam"
		betacontact_age[1:3,:]=[3.424 1.396 0.631 0.441 0.655 1.08 1.26 1.078 0.624 0.368 0.39 0.303 0.168 0.103 0.062 0.037 0.037 0.037 0.037 0.037;1.316 6.409 1.498 0.479 0.306 0.767 1.069 1.121 0.92 0.401 0.27 0.232 0.16 0.093 0.046 0.037 0.037 0.037 0.037 0.037;0.406 2.212 8.116 1.164 0.506 0.465 0.678 0.926 1.088 0.606 0.335 0.156 0.086 0.075 0.058 0.046 0.046 0.046 0.046 0.046];
		betacontact_age[4:6,:]=[0.258 0.608 3.35 11.788 2.213 1.134 0.848 1.088 1.242 1.139 0.614 0.253 0.1 0.051 0.03 0.022 0.022 0.022 0.022 0.022;0.414 0.345 0.485 3.478 5.118 2.538 1.595 1.409 1.168 1.284 0.908 0.471 0.147 0.038 0.041 0.035 0.035 0.035 0.035 0.035;0.856 0.492 0.329 1.261 2.945 4.328 2.415 1.858 1.52 1.23 1.165 0.616 0.223 0.054 0.024 0.022 0.022 0.022 0.022 0.022];
		betacontact_age[7:9,:]=[0.894 1.215 0.953 0.684 1.419 2.386 3.173 2.303 1.72 1.307 0.981 0.654 0.269 0.072 0.039 0.033 0.033 0.033 0.033 0.033;0.857 1.315 1.14 1.023 0.994 1.809 2.159 3.085 2.335 1.483 1.054 0.537 0.256 0.118 0.066 0.026 0.026 0.026 0.026 0.026;0.523 0.952 1.216 1.406 1.192 1.544 1.922 2.135 2.726 1.805 1.253 0.455 0.221 0.094 0.058 0.031 0.031 0.031 0.031 0.031];
		betacontact_age[10:12,:]=[0.462 0.693 0.856 1.766 1.116 1.256 1.459 1.652 1.697 1.892 1.177 0.543 0.164 0.056 0.052 0.053 0.053 0.053 0.053 0.053;0.467 0.767 1.072 1.487 1.282 1.659 1.466 1.394 1.8 1.864 1.663 0.879 0.24 0.064 0.052 0.055 0.055 0.055 0.055 0.055;0.772 0.894 0.783 0.992 0.888 1.42 1.464 1.129 1.225 1.026 1.156 0.992 0.356 0.114 0.05 0.044 0.044 0.044 0.044 0.044];
		betacontact_age[13:15,:]=[0.571 0.533 0.394 0.481 0.452 0.69 0.786 0.818 0.633 0.527 0.484 0.5 0.392 0.167 0.077 0.03 0.03 0.03 0.03 0.03;0.304 0.46 0.372 0.204 0.227 0.301 0.447 0.463 0.374 0.19 0.19 0.218 0.193 0.253 0.093 0.034 0.034 0.034 0.034 0.034;0.155 0.444 0.389 0.364 0.136 0.256 0.238 0.429 0.421 0.326 0.247 0.173 0.223 0.183 0.237 0.111 0.111 0.111 0.111 0.111];
		betacontact_age[16:end,:]=[0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131];
	end

		# Variable initial conditions

	S=zeros(nbr,Na);		# Susceptible
	Is=zeros(nbr,Na,Ni);	# Severe infected
	Im=zeros(nbr,Na,Ni);	# Mild infected
	Ip=zeros(nbr,Na,Ni);	# Paucisymptomatic infected
	H=zeros(nbr);			# Number of hospitalized at each time
	R=zeros(nbr,Na);		# Recovered population
	Mindir=zeros(nbr,Na);	# Population that died directly from COVID-19
	Mdir=zeros(nbr,Na);		# Population that died indirectly from COVID-19
	Ms=zeros(nbr,Na);		# Natural mortality (to compute proportion of infected in the end)
	NewHosp=zeros(nbr); 	# The inflow of new hospitalisations
	NewInfec=zeros(nbr);	# The inflow of new infections

	if country=="France"
		prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.03,0.02,0.01]; # age distribution of hospitalised (mean data between reports of Santé Publique France and CDC
		Is[1,:,1:NiS]=(130/(NiS*di))*prop_age*ones(NiS)'; 	# 130 is the number of severe infected in France at March, 1st, uniformly distributed in age since infection then distributed in age
		Im[1,:,1:NiS]=(130/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
		Ip[1,:,1:NiS]=(130/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
		S[1,:]=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
	elseif country=="Burkina"
		prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.06]; #
		Is[1,:,1:NiS]=(288/(NiS*di))*prop_age*ones(NiS)'; 	# 288 is the number of cases in Burkina, on April, 1st (WHO ref)
		Im[1,:,1:NiS]=(288/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
		Ip[1,:,1:NiS]=(288/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
		S[1,:]=18450494/100*[13.2,12.6,11.6,12.9,11.5,9.3,7.2,5.4,4.3,3.1,2.5,1.8,1.4,0.9,0.7,0.4,0.3,0.2]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
	elseif country=="Vietnam"
		prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.03,0.02,0.01]; # age distribution of hospitalised (mean data between reports of Santé Publique France and CDC
		Is[1,:,1:NiS]=(217/(NiS*di))*prop_age*ones(NiS)'; 	# 130 is the number of severe infected in France at March, 1st, uniformly distributed in age since infection then distributed in age
		Im[1,:,1:NiS]=(217/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
		Ip[1,:,1:NiS]=(217/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
		S[1,:]=[4539031,7525542,6976064,6474991,7145151,8741274,8345522,7631483,6941367,6442420,5774836,5136526,4149565,2772903,1526235,1113819,873832,606471,268159,124992]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
	end

	H[1]=di*sum(Is[1,:,Tsympt:end]);

	muadd=(0.01*munat)/(1+99*exp(-10*(((H[1])/(Hsat))-1)));
	gammaIindir=(gammaIdir)/(1+99*exp(-10*(((H[1])/(Hsat))-1)));

	Mdir[1,:]=sum(gammaIdir.*Is[1,:,:],dims=2)[:,1];
	Mindir[1,:]=sum(gammaIindir.*Is[1,:,:],dims=2)[:,1]+muadd.*S[1,:]+di*sum(muadd.*Is[1,:,:],dims=2)[:,1]+di*sum(muadd.*Im[1,:,:],dims=2)[:,1]+di*sum(muadd.*Ip[1,:,:],dims=2)[:,1]+muadd.*R[1,:];
	Ms[1,:]=muadd.*S[1,:]+munat.*S[1,:];
	NewHosp[1]=di*sum(Is[1,:,:]);
	NewInfec[1]=di*sum(Is[1,:,:])+di*sum(Im[1,:,:])+di*sum(Ip[1,:,:]);

		# R0 computation

 	R0=0;
	# piS=zeros(Na,Ni);
	# piM=zeros(Na,Ni);
	# pIp=zeros(Na,Ni);
	# omega=zeros(Na,Ni);
	#
	# for aa=1:Na
	# 	piS[aa,:]=exp.(-munat[aa]*i-(i-imaxS*ones(Ni)).*Charac.(i,imaxS,imax+di)-gammaIdir[aa,:].*(i-isympt*ones(Ni)).*Charac.(i,isympt,imaxS));
	# 	piM[aa,:]=exp.(-munat[aa]*i-(i-imaxM*ones(Ni)).*Charac.(i,imaxM,imax+di));
	# 	pIp[aa,:]=exp.(-munat[aa]*i-(i-imaxM*ones(Ni)).*Charac.(i,imaxM,imax+di));
	# 	omega[aa,:]=(1-p)*prophosp[aa]*betaS.*piS[aa,:]+(1-p)*(1-prophosp[aa])*betaM.*piM[aa,:]+p*betaA.*pIp[aa,:];
	# end
	#
	# Omega=zeros(Na)
	# Omega=di*sum(omega,dims=2)[:,1];
	# Omegamat=hcat(fill.(Omega,Na)...)';
	#
	# S0mat=hcat(fill.(S[1,:],Na)...)';
	# R0mat=betacontact_age.*Omegamat.*S0mat;
	# R0=maximum(abs.(eigvals(R0mat)));		# need package LinearAlgebra to compute eigenvalues

		# Confinement

	betacontact=ones(Na,Na);			# contract probability updated (according to control measure).

		### Boucle ###

	S=convert(Array{Float64,2},S); 	# use typeof(S) to know its type, so that we need to convert to have no problem
	Is=convert(Array{Float64,3},Is);
	Im=convert(Array{Float64,3},Im);
	Ip=convert(Array{Float64,3},Ip);

	Lambda=zeros(nbr,Na);			# force of infection (age dependent)

	for tt=1:(nbr-1)

		betacontact=betacontact_age.*(ones(Na,Na)-hcat(fill.(c[tt,:],Na)...)'); # repeat 10 times the column vector c[tt,:] horizontally, which we multiply by vector of contacts by age

			#  We compute the mortality rates according to the number of infected (according to the the saturation threshold Hsat of hospitals)

		muadd=(0.01*munat)/(1+99*exp.(-10*(((H[tt])/(Hsat))-1)));
		gammaIindir=(gammaIdir)/(1+99*exp(-10*(((H[tt])/(Hsat))-1)));

			# On continue la boucle ##

		Lambda1=di*Is[tt,:,:]*betaS+di*Im[tt,:,:]*betaM+di*Ip[tt,:,:]*betaA; 	# we compute the force of infection for each class of age
		Lambda[tt,:]=betacontact*Lambda1;

		S[tt+1,:]=(S[tt,:])./(ones(Na)+dt*munat+dt*muadd+dt*Lambda[tt,:]); 		# we compute the new density of susceptible
		Is[tt+1,:,1]=(Is[tt,:,1]+(1-p)*(dt/di)*prophosp.*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hS[:,1]+(dt/di)*gammaIdir[:,1]+(dt/di)*gammaIindir[:,1]);	# we compute the new density of infected that have been contaminated
		Im[tt+1,:,1]=(Im[tt,:,1]+(1-p)*(dt/di)*(ones(Na)-prophosp).*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,1]); # mild
		Ip[tt+1,:,1]=(Im[tt,:,1]+p*(dt/di)*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,1]); 		# paucisymptomatic

		for ii=1:(Ni-1)
			Is[tt+1,:,ii+1]=(Is[tt,:,ii+1]+(dt/di)*Is[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hS[:,ii+1]+(dt/di)*gammaIdir[:,ii+1]+(dt/di)*gammaIindir[:,ii+1]);
			Im[tt+1,:,ii+1]=(Im[tt,:,ii+1]+(dt/di)*Im[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,ii+1]);
			Ip[tt+1,:,ii+1]=(Ip[tt,:,ii+1]+(dt/di)*Ip[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,ii+1]);
		end

		H[tt+1]=di*sum(Is[tt+1,:,Tsympt:end]);

		NewHosp[tt+1]=sum((1-p)*prophosp.*Lambda[tt,:].*S[tt+1,:]);
		NewInfec[tt+1]=sum(Lambda[tt,:].*S[tt+1,:]);

		Mdir[tt+1,:]=Mdir[tt,:]+dt*sum(gammaIdir.*Is[tt+1,:,:],dims=2)[:,1];
		Mindir[tt+1,:]=Mindir[tt,:]+dt*sum(gammaIindir.*Is[tt+1,:,:],dims=2)[:,1]+dt*muadd.*S[tt+1,:]+dt*di*sum(muadd.*Is[tt+1,:,:],dims=2)[:,1]+dt*di*sum(muadd.*Im[tt+1,:,:],dims=2)[:,1]+dt*di*sum(muadd.*Ip[tt+1,:,:],dims=2)[:,1]+dt*muadd.*R[tt+1,:];
		R[tt+1,:]=(R[tt,:]+dt*Is[tt,:,Ni]+dt*Im[tt,:,Ni]+dt*Ip[tt,:,Ni]+dt*di*sum(hS.*Is[tt+1,:,:],dims=2)[:,1]+dt*di*sum(hM.*Im[tt+1,:,:],dims=2)[:,1]+dt*di*sum(hM.*Ip[tt+1,:,:],dims=2)[:,1])./(ones(Na)+dt*munat+dt*muadd);				# we compute the new number of individuals that have recovered from COVID-19
		Ms[tt,:]=dt*muadd.*S[tt,:]+dt*munat.*S[tt,:];

	end

	Hosp=zeros(nbr,Na);
	NonHosp=zeros(nbr,Na);

	Hosp=di*sum(Is[:,:,Tsympt:end],dims=3)[:,:,1];
	NonHosp=di*sum(Is[:,:,1:(Tsympt-1)],dims=3)[:,:,1]+di*sum(Im[:,:,:],dims=3)[:,:,1]+di*sum(Ip[:,:,:],dims=3)[:,:,1];

	prop_age=zeros(nbr,Na);		# proportion, by age, of infected, in time
	prop_totale=zeros(nbr);		# mean proportion of infected along time

	prop_age[1,:]=ones(Na)-(S[1,:]./(S[1,:]+di*sum(Is,dims=3)[1,:,1]+di*sum(Im,dims=3)[1,:,1]+di*sum(Ip,dims=3)[1,:,1]));
	prop_totale[1]=1-(sum(S[1,:])/(sum(S[1,:])+di*sum(sum(Is,dims=2),dims=3)[1,1,1]+di*sum(sum(Im,dims=2),dims=3)[1,1,1]+di*sum(sum(Ip,dims=2),dims=3)[1,1,1]));

	for tt=2:nbr
		prop_age[tt,:]=ones(Na)-(((S[tt,:])+sum(Ms[1:(tt-1),:],dims=1)[1,:])./(S[1,:]+di*sum(Is,dims=3)[1,:,1]+di*sum(Im,dims=3)[1,:,1]+di*sum(Ip,dims=3)[1,:,1]));
		prop_totale[tt]=1-((sum(S[tt,:])+sum(Ms[1:(tt-1)]))/(sum(S[1,:])+di*sum(sum(Is,dims=2),dims=3)[1,1,1]+di*sum(sum(Im,dims=2),dims=3)[1,1,1]+di*sum(sum(Ip,dims=2),dims=3)[1,1,1]));
	end

	## Figures of mortality rates

	if country=="France"
		HH=collect(0:100:20000);
	elseif country=="Burkina"
		HH=collect(0:1:100);
	elseif country=="Vietnam"
		HH=collect(0:100:20000);
	end

	nbrHH=length(HH);
	muHH=zeros(nbrHH,Na);
	gammaHH=zeros(nbrHH,Na);

 gammamort_age=[0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.025,0.0295,0.0335,0.0511,0.0711,0.117,0.157,0.196,0.236,0.2945,0.3245,0.3545,0.3845];

	for aa=1:Na
		muHH[:,aa]=munat[aa]*ones(nbrHH)+(0.01*munat[aa])./(ones(nbrHH)+99*exp.(-10*((HH/Hsat)-ones(nbrHH))));
		gammaHH[:,aa]=gammamort_age[aa]*ones(nbrHH)+(gammamort_age[aa])./(ones(nbrHH)+99*exp.(-10*((HH/Hsat)-ones(nbrHH))));
	end

	pltmu=plot(HH,a*5,muHH',st=:surface,title="Death rate \$\\mu\$(a,H)",camera=(-60,30),yticks=collect(0:10:100),xtickfont=font(26),ytickfont=font(26),ztickfont=font(26),xguidefontsize=30,yguidefontsize=30,zguidefontsize=30,titlefontsize=34,colorbar=false);
	pltgamma=plot(HH,a*5,gammaHH',st=:surface,title="Mortality rate \$\\gamma\$(a,H)",camera=(-60,30),yticks=collect(0:10:100),xtickfont=font(26),ytickfont=font(26),ztickfont=font(26),xguidefontsize=30,yguidefontsize=30,zguidefontsize=30,titlefontsize=34,colorbar=false);

	return t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma;
end

	##### Epidemic scenario without control in each country (France, Burkina and Vietnam) - Figure 4 #####

function CovidFig(T,country)

	timecomput=time();	# initial time

	if country=="France"
		Hsat=5000;
		amax=20;
	elseif country=="Burkina"
		Hsat=11;
		amax=18;
	elseif country=="Vietnam"
		amax=20;
		Hsat=5932;
	else
		println("You did not write the country correctly, please type either France, Burkina, or Vietnam")
		return
	end

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma=Covid(T,0,Hsat,0.5,5.2,0.1,country);

	Na=length(a);
	tt=ones(length(t));

	## Toutes les figures de S,H,NH,R,M

	pp=plot(t,log10.(tt+sum(S,dims=2)[:,1]),title="Populations over time",color="maroon",xlabel="Time (days)",xtickfont=font(26),ytickfont=font(26),xguidefontsize=30,yguidefontsize=30,ylabel="\$\\log_{10}\$(populations)",label="Susceptible",titlefontsize=34,width=2);
	pp=plot!(pp,t,log10.(tt+sum(Hosp,dims=2)[:,1]),label="Hospitalised",color="orange",width=2);
	pp=plot!(pp,t,log10.(tt+sum(NonHosp,dims=2)[:,1]),label="Non hospitalised",color="darkgreen",width=2);
	pp=plot!(pp,t,log10.(tt+sum(R,dims=2)[:,1]),label="Recovered",color="darkblue",width=2);
	pp=plot!(pp,t,log10.(tt+sum(Mdir+Mindir,dims=2)[:,1]),label="Deceased",color="red",width=2);
	pp=plot!(pp,t,log10.(tt+Hsat*ones(length(t))),linestyle=:dash,label="Healthcare capacity",color="orange",width=2);
	pp=plot!(pp,xticks = (collect(0:15:T)),legendfont = font("Times new roman", 24),yticks=collect(0:1:8),legend=:bottom,width=2);
	pp=plot!(pp,left_margin=-1.1cm, right_margin=-0.1cm,bottom_margin=-0.9cm,top_margin=-0.3cm);

	## Figures en 3D, de Morts, Hospitalisés et NonHospitalisés

	pltMage=plot(t,a*5,(Mdir+Mindir)',st=:surface,camera=(-60,30),title="Number of cumulative deaths",xtickfont=font(26),ytickfont=font(26), ztickfont=font(26),xguidefontsize=30,yguidefontsize=30,zguidefontsize=30,legend=false,titlefontsize=34);
	pltMage=plot!(pltMage,xticks=collect(0:25:T),yticks=collect(0:10:100),color="red",legendfont = font("Times new roman", 24));

	pltH=plot(t,a*5,Hosp',st=:surface,camera=(-60,30),title="Number of hospitalised", xtickfont=font(26),ytickfont=font(26),ztickfont=font(26),xguidefontsize=30,yguidefontsize=30,zguidefontsize=30,legend=false,titlefontsize=34);
	pltH=plot!(pltH,xticks=collect(0:25:T),yticks=collect(0:10:100),legendfont = font("Times new roman", 24));

	pltNH=plot(t,a*5,NonHosp',st=:surface,camera=(-60,30),title="Number of non hospitalised",xtickfont=font(26),ytickfont=font(26),ztickfont=font(26),xguidefontsize=30,yguidefontsize=30,zguidefontsize=30,legend=false,titlefontsize=34);
	pltNH=plot!(pltNH,xticks=collect(0:25:T),yticks=collect(0:10:100),legendfont = font("Times new roman", 24));

	## Figures de la proportion d'infectés selon la classe d'âge

	proptot=prop_totale[end];
	propage=prop_age[end,:];
	plotpropage=plot(collect(0:5:(amax)*5),[propage;propage[end]],xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xlabel="Age (years)",ylabel="Proportion",title="Proportion of infected at t=$T",label="Proportion by age class",titlefontsize=38,ylims=(0,1),width=2,linetype=:steppost);
	plotpropage=plot!(plotpropage,xticks = (collect(0:5:(amax*5))),yticks=collect(0:0.1:1),legendfont = font("Times new roman", 28));
	plotpropage=plot!(plotpropage,collect(0:0.1:(amax)*5),proptot*ones(length(collect(0:0.1:((amax-1)*5)))),linestyle = :dash,label="Average proportion",color="red",width=2);
	plotpropage=plot!(plotpropage,collect(0:0.1:((amax)*5)),(1-(1/3.3))*ones(Na),linestyle = :dash,label="Herd immunity",color="green",legend=:bottomright,width=2);
	plotpropage=plot!(plotpropage,left_margin=-1.9cm, right_margin=-2.5cm,bottom_margin=-0.9cm,top_margin=-0.3cm);

	println("The function took $(time()-timecomput) seconds to execute.")

	return pp,pltMage,pltH,pltNH,plotpropage;
end

	##### Age-structure of each country - Figure 3 - and Proportion of infected at the end of the epidemics without measure controls - Figure 5 -   #####

function ComparisonCountries(T)

	amin=1;amax=20;			# discretisation in age
	Na=Int(amax); 			# number of groups (20)
	a=collect(amin:1:amax); # vector of age

		# Age-structure populations

	PopFrance=zeros(Na);
	PopBurkina=zeros(Na);
	PopVietnam=zeros(Na);

	PopFrance=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655];
	PopBurkina=18450494/100*[13.2,12.6,11.6,12.9,11.5,9.3,7.2,5.4,4.3,3.1,2.5,1.8,1.4,0.9,0.7,0.4,0.3,0.2];
	PopVietnam=[4539031,7525542,6976064,6474991,7145151,8741274,8345522,7631483,6941367,6442420,5774836,5136526,4149565,2772903,1526235,1113819,873832,606471,268159,124992];

	pltPop=plot(collect(0:5:(amax-2)*5),[PopBurkina;PopBurkina[end]],xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xlabel="Age (years)",ylabel="Population",title="Age-structure of the populations",label="Burkina Faso",titlefontsize=38,width=2,linetype=:steppost,color="blue");
	pltPop=plot!(pltPop,xticks = (collect(0:5:(amax*5))),legendfont = font("Times new roman", 28));
	pltPop=plot!(pltPop,left_margin=-1.9cm, right_margin=-2.5cm,bottom_margin=-0.9cm,top_margin=-0.3cm);
	pltPop=plot!(pltPop,collect(0:5:(amax)*5),[PopFrance;PopFrance[end]],label="France",titlefontsize=38,width=2,linetype=:steppost,color="limegreen");
	pltPop=plot!(pltPop,collect(0:5:(amax)*5),[PopVietnam;PopVietnam[end]],label="Vietnam",titlefontsize=38,width=2,linetype=:steppost,color="darkorange");

	pltPopProp=plot(collect(0:5:(amax-2)*5),[PopBurkina/sum(PopBurkina);PopBurkina[end]/sum(PopBurkina)],xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xlabel="Age (years)",ylabel="Frequency",title="Age-frequency of the populations",label="Burkina Faso",titlefontsize=38,width=2,linetype=:steppost,color="blue");
	pltPopProp=plot!(pltPopProp,xticks = (collect(0:5:(amax*5))),yticks=collect(0:0.02:1),legendfont = font("Times new roman", 28));
	pltPopProp=plot!(pltPopProp,left_margin=-1.9cm, right_margin=-2.5cm,bottom_margin=-0.9cm,top_margin=-0.3cm);
	pltPopProp=plot!(pltPopProp,collect(0:5:(amax)*5),[PopFrance/sum(PopFrance);PopFrance[end]/sum(PopFrance)],label="France",titlefontsize=38,width=2,linetype=:steppost,color="limegreen");
	pltPopProp=plot!(pltPopProp,collect(0:5:(amax)*5),[PopVietnam/sum(PopVietnam);PopVietnam[end]/sum(PopVietnam)],label="Vietnam",titlefontsize=38,width=2,linetype=:steppost,color="darkorange");

		 # Age distribution of the proportion at the end of the epidemics (without control measures)

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma=Covid(T,0,11,0.5,5.2,0.1,"Burkina");
	propage=prop_age[end,:];
	proptot=prop_totale[end];

	plotprop=plot(collect(0:5:(amax-2)*5),[propage;propage[end]],xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xlabel="Age (years)",ylabel="Proportion",title="Proportion of infected at t=$T",label="Burkina Faso",titlefontsize=38,ylims=(0,1),width=2,linetype=:steppost,color="blue");
	plotprop=plot!(plotprop,xticks = (collect(0:5:(amax*5))),yticks=collect(0:0.1:1),legendfont = font("Times new roman", 28));
	plotprop=plot!(plotprop,collect(0:0.1:(amax-2)*5),proptot*ones(length(collect(0:0.1:((amax-3)*5)))),linestyle = :dash,label="",color="blue",width=2);

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma=Covid(T,0,5000,0.5,5.2,0.1,"France");
	propage=prop_age[end,:];
	proptot=prop_totale[end];

	plotprop=plot!(plotprop,collect(0:5:(amax)*5),[propage;propage[end]],label="France",width=2,linetype=:steppost,color="limegreen");
	plotprop=plot!(plotprop,collect(0:0.1:(amax)*5),proptot*ones(length(collect(0:0.1:((amax-1)*5)))),linestyle = :dash,label="",color="limegreen",width=2);

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma=Covid(T,0,5932,0.5,5.2,0.1,"Vietnam");
	propage=prop_age[end,:];
	proptot=prop_totale[end];

	plotprop=plot!(plotprop,collect(0:5:(amax)*5),[propage;propage[end]],label="Vietnam",width=2,linetype=:steppost,color="darkorange");
	plotprop=plot!(plotprop,collect(0:0.1:(amax)*5),proptot*ones(length(collect(0:0.1:((amax-1)*5)))),linestyle = :dash,label="",color="darkorange",width=2);

	plotprop=plot!(plotprop,collect(0:0.1:((amax)*5)),(1-(1/3.3))*ones(Na),linestyle = :dash,label="Herd immunity",color="black",legend=:bottomright,width=2);
	plotprop=plot!(plotprop,left_margin=-1.9cm, right_margin=-2.5cm,bottom_margin=-0.9cm,top_margin=-0.3cm);
	plotprop=plot!(plotprop,legend=:bottom)

		# Social contacts matrices

	## Data from Beraud, Boelle et al (2015)

	# amax=15;
	# ContactFr=zeros(amax,amax);
	# ContactFr[1:3,:]=[4.42e-07 3.40e-07 1.21e-07 5.53e-08 4.39e-08 9.27e-08 2.06e-07 1.75e-07 1.10e-07 8.05e-08 5.01e-08 6.02e-08 6.76e-08 4.43e-08 2.52e-08; 3.40e-07 8.25e-07 5.71e-07 1.58e-07 5.03e-08 5.67e-08 1.11e-07 1.75e-07 2.11e-07 1.62e-07 6.37e-08 4.60e-08 5.56e-08 6.85e-08 3.43e-08;1.21e-07 5.71e-07 1.63e-06 9.87e-07 1.46e-07 5.63e-08 5.65e-08 9.46e-08 2.29e-07 1.96e-07 9.69e-08 5.17e-08 3.53e-08 4.94e-08 4.80e-08];
	# ContactFr[4:6,:]=[5.53e-08 1.58e-07 9.87e-07 1.28e-06 8.56e-07 2.35e-07 8.02e-08 7.58e-08 1.35e-07 2.11e-07 1.82e-07 1.08e-07 5.50e-08 4.44e-08 3.71e-08; 4.39e-08 5.03e-08 1.46e-07 8.56e-07 7.13e-07 4.96e-07 2.23e-07 1.49e-07 1.22e-07 1.76e-07 1.73e-07 1.64e-07 1.17e-07 4.99e-08 4.58e-08; 9.27e-08 5.67e-08 5.63e-08 2.35e-07 4.96e-07 5.16e-07 4.22e-07 2.63e-07 1.79e-07 2.02e-07 1.76e-07 2.08e-07 1.65e-07 7.26e-08 5.79e-08];
	# ContactFr[7:9,:]=[2.06e-07 1.11e-07 5.65e-08 8.02e-08 2.23e-07 4.22e-07 3.46e-07 3.43e-07 2.45e-07 2.38e-07 1.98e-07 1.94e-07 2.11e-07 1.30e-07 7.38e-08; 1.75e-07 1.75e-07 9.46e-08 7.58e-08 1.49e-07 2.63e-07 3.43e-07 3.70e-07 3.26e-07 2.32e-07 1.77e-07 1.59e-07 1.56e-07 1.45e-07 1.35e-07; 1.10e-07 2.11e-07 2.29e-07 1.35e-07 1.22e-07 1.79e-07 2.45e-07 3.26e-07 3.11e-07 2.94e-07 2.31e-07 1.86e-07 1.45e-07 1.22e-07 9.07e-08];
	# ContactFr[10:12,:]=[8.05e-08 1.62e-07 1.96e-07 2.11e-07 1.76e-07 2.02e-07 2.38e-07 2.32e-07 2.94e-07 3.24e-07 3.05e-07 2.53e-07 1.97e-07 1.33e-07 7.50e-08; 5.01e-08 6.37e-08 9.69e-08 1.82e-07 1.73e-07 1.76e-07 1.98e-07 1.77e-07 2.31e-07 3.05e-07 3.28e-07 3.13e-07 2.15e-07 1.19e-07 8.92e-08; 6.02e-08 4.60e-08 5.17e-08 1.08e-07 1.64e-07 2.08e-07 1.94e-07 1.59e-07 1.86e-07 2.53e-07 3.13e-07 3.29e-07 2.52e-07 1.45e-07 9.38e-08];
	# ContactFr[13:15,:]=[6.76e-08 5.56e-08 3.53e-08 5.50e-08 1.17e-07 1.65e-07 2.11e-07 1.56e-07 1.45e-07 1.97e-07 2.15e-07 2.52e-07 2.36e-07 2.28e-07 1.69e-07; 4.43e-08 6.85e-08 4.94e-08 4.44e-08 4.99e-08 7.26e-08 1.30e-07 1.45e-07 1.22e-07 1.33e-07 1.19e-07 1.45e-07 2.28e-07 3.60e-07 2.49e-07; 2.52e-08 3.43e-08 4.80e-08 3.71e-08 4.58e-08 5.79e-08 7.38e-08 1.35e-07 9.07e-08 7.50e-08 8.92e-08 9.38e-08 1.69e-07 2.49e-07 2.13e-07];

	## Data from Prem et al (2017)

	amax=16;
	ContactFr=zeros(amax,amax);
	ContactFr=[3.805 1.106 0.412 0.255 0.342 0.753 1.209 1.042	 0.496	 0.31	 0.307	 0.272	 0.186	 0.128	 0.088	 0.054;
	1.062	 5.033	 1.011	 0.273	 0.164	 0.415	 0.921	 1.138	 0.845	 0.331	 0.207	 0.185	 0.177	 0.112	 0.058	 0.053;
	0.238	 1.538	 6.986	 0.846	 0.287	 0.244	 0.403	 0.827	 1.081	 0.57	 0.266	 0.147	 0.102	 0.088	 0.07	 0.07;
	0.124	 0.308	 2.301	 7.832	 1.36	 0.651	 0.531	 0.762	 1.025	 1.049	 0.479	 0.164	 0.072	 0.053	 0.035	 0.026;
	0.2	 0.168	 0.256	 2.165	 3.934	 1.708	 1.168	 1.133	 0.97	 1.27	 0.813	 0.327	 0.084	 0.046	 0.054	 0.052;
	0.533	 0.249	 0.144	 0.731	 1.944	 3.458	 1.777	 1.471	 1.263	 0.997	 0.961	 0.357	 0.113	 0.052	 0.031	 0.024;
	0.722	 0.848	 0.53	 0.423	 1.005	 1.72	 2.919	 1.894	 1.457	 1.135	 0.81	 0.405	 0.184	 0.085	 0.047	 0.053;
	0.704	 1.072	 0.837	 0.683	 0.765	 1.407	 1.76	 3.211	 2.098	 1.344	 0.92	 0.364	 0.26	 0.152	 0.098	 0.046;
	0.31	 0.657	 0.983	 1.109	 0.924	 1.252	 1.591	 1.877	 2.954	 1.682	 1.113	 0.273	 0.184	 0.111	 0.087	 0.048;
	0.406	 0.467	 0.612	 1.537	 0.945	 0.994	 1.259	 1.465	 1.573	 2.196	 1.086	 0.341	 0.151	 0.086	 0.083	 0.082;
	0.252	 0.595	 0.829	 1.202	 1.019	 1.291	 1.179	 1.215	 1.594	 1.743	 1.918	 0.649	 0.267	 0.119	 0.088	 0.09;
	0.542	 0.671	 0.574	 0.749	 0.614	 0.932	 0.918	 0.748	 0.793	 0.681	 0.923	 1.462	 0.513	 0.218	 0.12	 0.095;
	0.396	 0.363	 0.257	 0.361	 0.339	 0.508	 0.639	 0.694	 0.523	 0.416	 0.441	 0.675	 1.406	 0.426	 0.263	 0.129;
	0.198	 0.313	 0.257	 0.157	 0.212	 0.307	 0.45	 0.454	 0.425	 0.302	 0.324	 0.44	 0.518	 1.003	 0.276	 0.139;
	0.105	 0.294	 0.308	 0.355	 0.156	 0.245	 0.253	 0.477	 0.515	 0.434	 0.339	 0.338	 0.66	 0.619	 0.974	 0.326;
	0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643];


	my_cg = cgrad([:black,:darkblue,:blue,:darkgreen,:green,:yellow,:orange,:darkorange,:red]);
	pltFr=plot(collect(0:5:(amax-1)*5),collect(0:5:(amax-1)*5),c=my_cg,log.(ContactFr),st=:heatmap,xtickfont=font(30),ytickfont=font(30),xticks=collect(0:10:amax*5),yticks=collect(0:10:amax*5),xguidefontsize=34,yguidefontsize=34,xlabel="Age of individuals (years)",ylabel="Age of contact",title="France",titlefontsize=38);

	amax=16;
	ContactBk=zeros(16,16);
	ContactBk[1:3,:]=[6.415 3.639 1.748 0.919 1.198 1.574 1.469 1.158 0.661 0.398 0.433 0.377 0.251 0.163 0.079 0.036;3.477 18.279 4.124 1.077 0.652 1.292 1.3506 1.237 0.911 0.452 0.356 0.345 0.283 0.160 0.068 0.037;1.177 5.907 17.935 1.868 0.832 0.717 0.830 0.862 0.851 0.500 0.314 0.169 0.1295 0.113 0.068 0.033];
	ContactBk[4:6,:]=[0.634 1.532 5.808 13.816 2.361 1.144 0.784 0.923 0.852 0.703 0.440 0.228 0.141 0.064 0.031 0.015;0.849 0.858 0.911 4.087 4.586 2.153 1.273 1.049 0.824 0.768 0.581 0.346 0.179 0.038 0.029 0.017;1.172 0.782 0.525 1.225 2.435 2.846 1.660 1.282 1.038 0.797 0.725 0.407 0.221 0.038 0.013 0.010];
	ContactBk[7:9,:]=[1.125 1.647 1.311 0.664 1.163 1.753 1.842 1.502 1.175 0.902 0.657 0.471 0.248 0.043 0.021 0.013;0.900 1.390 1.052 0.789 0.725 1.249 1.356 1.675 1.461 0.962 0.744 0.385 0.200 0.063 0.028 0.009;0.738 1.218 1.200 1.010 0.833 1.094 1.302 1.365 1.516 1.168 0.891 0.366 0.255 0.077 0.029 0.013];
	ContactBk[10:12,:]=[0.464 0.922 0.951 1.188 0.690 0.873 1.009 1.117 1.070 1.035 0.781 0.448 0.186 0.044 0.027 0.022;0.565 0.920 1.185 1.005 0.724 1.029 1.013 0.980 1.238 1.241 1.001 0.639 0.234 0.040 0.030 0.022;1.290 1.653 1.259 0.956 0.700 1.103 1.183 0.908 0.985 0.923 0.928 0.627 0.346 0.108 0.033 0.038];
	ContactBk[13:15,:]=[0.827 0.941 0.750 0.525 0.443 0.629 0.684 0.775 0.693 0.622 0.575 0.506 0.219 0.098 0.039 0.010;0.575 0.935 0.720 0.257 0.232 0.246 0.335 0.357 0.284 0.163 0.189 0.174 0.126 0.083 0.059 0.017;0.244 0.792 0.6343 0.393 0.103 0.205 0.141 0.258 0.241 0.238 0.243 0.144 0.154 0.099 0.061 0.061];
	ContactBk[16,:]=[0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022];

	pltBk=plot(collect(0:5:(amax-1)*5),collect(0:5:(amax-1)*5),c=my_cg,log.(ContactBk),st=:heatmap,xtickfont=font(30),ytickfont=font(30),xticks=collect(0:10:amax*5),yticks=collect(0:10:amax*5),xguidefontsize=34,yguidefontsize=34,xlabel="Age of individuals (years)",ylabel="Age of contact",title="Burkina Faso",titlefontsize=38,legend=true);

	amax=16;

	ContactViet=zeros(16,16);
	ContactViet[1:3,:]=[3.424 1.396 0.631 0.441 0.655 1.08 1.26 1.078 0.624 0.368 0.39 0.303 0.168 0.103 0.062 0.037;1.316 6.409 1.498 0.479 0.306 0.767 1.069 1.121 0.92 0.401 0.27 0.232 0.16 0.093 0.046 0.037;0.406 2.212 8.116 1.164 0.506 0.465 0.678 0.926 1.088 0.606 0.335 0.156 0.086 0.075 0.058 0.046];
	ContactViet[4:6,:]=[0.258 0.608 3.35 11.788 2.213 1.134 0.848 1.088 1.242 1.139 0.614 0.253 0.1 0.051 0.03 0.022;0.414 0.345 0.485 3.478 5.118 2.538 1.595 1.409 1.168 1.284 0.908 0.471 0.147 0.038 0.041 0.035;0.856 0.492 0.329 1.261 2.945 4.328 2.415 1.858 1.52 1.23 1.165 0.616 0.223 0.054 0.024 0.022];
	ContactViet[7:9,:]=[0.894 1.215 0.953 0.684 1.419 2.386 3.173 2.303 1.72 1.307 0.981 0.654 0.269 0.072 0.039 0.033;0.857 1.315 1.14 1.023 0.994 1.809 2.159 3.085 2.335 1.483 1.054 0.537 0.256 0.118 0.066 0.026;0.523 0.952 1.216 1.406 1.192 1.544 1.922 2.135 2.726 1.805 1.253 0.455 0.221 0.094 0.058 0.031];
	ContactViet[10:12,:]=[0.462 0.693 0.856 1.766 1.116 1.256 1.459 1.652 1.697 1.892 1.177 0.543 0.164 0.056 0.052 0.053;0.467 0.767 1.072 1.487 1.282 1.659 1.466 1.394 1.8 1.864 1.663 0.879 0.24 0.064 0.052 0.055;0.772 0.894 0.783 0.992 0.888 1.42 1.464 1.129 1.225 1.026 1.156 0.992 0.356 0.114 0.05 0.044];
	ContactViet[13:15,:]=[0.571 0.533 0.394 0.481 0.452 0.69 0.786 0.818 0.633 0.527 0.484 0.5 0.392 0.167 0.077 0.03;0.304 0.46 0.372 0.204 0.227 0.301 0.447 0.463 0.374 0.19 0.19 0.218 0.193 0.253 0.093 0.034;0.155 0.444 0.389 0.364 0.136 0.256 0.238 0.429 0.421 0.326 0.247 0.173 0.223 0.183 0.237 0.111];
	ContactViet[16,:]=[0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131];

	pltViet=plot(collect(0:5:(amax-1)*5),collect(0:5:(amax-1)*5),c=my_cg,log.(ContactViet),st=:heatmap,xtickfont=font(30),ytickfont=font(30),xticks=collect(0:10:amax*5),yticks=collect(0:10:amax*5),xguidefontsize=34,yguidefontsize=34,xlabel="Age of individuals (years)",ylabel="Age of contact",title="Vietnam",titlefontsize=38);

	return pltPop,pltPopProp,plotprop,pltFr,pltBk,pltViet;
end

	#####  Sensibility analysis, data for Figure S1 #####

function SensAnal(T)

	timecomput=time();	# initial time

	Hsat=[10,100,500,2000,5000,50000,500000,5000000];
	p=collect(0.05:0.1:0.95);
	isympt=collect(1.2:2:9.2);
	xi_p=[0.1,0.3,0.5,0.7,1];
	country=["Burkina","France","Vietnam"];
	size_hsat=length(Hsat);
	size_p=length(p);
	size_isympt=length(isympt);
	size_xip=length(xi_p);
	size_country=length(country)
	sizetotal=size_hsat*size_p*size_isympt*size_xip*size_country;

	SensAnalMatrix=zeros(sizetotal,8);
	zz=1;

	for gg=1:size_country
		for hh=1:size_xip
			for ii=1:size_hsat
				for jj=1:size_p
					for kk=1:size_isympt
						t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale,R0,pltmu,pltgamma=Covid(T,0,Hsat[ii],p[jj],isympt[kk],xi_p[hh],country[gg]);
						SensAnalMatrix[zz,1]=sum(Mdir[end,:]+Mindir[end,:]);
						SensAnalMatrix[zz,2]=0.2*sum(NewHosp);
						SensAnalMatrix[zz,3]=0.2*sum(NewInfec);
						SensAnalMatrix[zz,4]=gg;
						SensAnalMatrix[zz,5]=xi_p[hh];
						SensAnalMatrix[zz,6]=Hsat[ii];
						SensAnalMatrix[zz,7]=p[jj];
						SensAnalMatrix[zz,8]=isympt[kk];
						zz=zz+1;
					end
				end
			end
		end
	end

	println("Deaths|Hosp|Infec|Country|xip|Hsat|p|isympt")
	println("The function took $(time()-timecomput) seconds to execute.")

	return SensAnalMatrix;
end

############################# OPTIMAL CONTROL ##############################

		##### PDE model + Adjoint system #####

function Covidadj(T,c,B,p,country)

		# Time;

	dt=0.2;				# discretisation time (in days)
	nbr=Int(T/dt)+1;	# number of iterations
	t=collect(0:dt:T);	# vector of time

		# Age

	amin=1;		# discretisation in age

	if country=="France"
		amax=20;	# (20 groups of 5 years age, first one is 0-5 years, 2nd: 5-10, etc. - last one is 95 and more)
		Hsat=5000;	# 5000 ICU beds in France
	elseif country=="Burkina"
		amax=18;	# (18 groups of 5 year age, - last one is 85 and more-)
		Hsat=11; 	# 11 in Burkina Faso
	elseif country=="Vietnam"
		amax=20;
		Hsat=5932;
	else
		println("You did not write the country correctly, type either France, Burkina, or Vietnam")
		return
	end


	Na=Int(amax);
	a=collect(amin:1:amax);

	if length(c)==1
		c=c*ones(nbr,Na);
	end
	if length(B)==1
		B=B*ones(Na);
	end

		# Age of infection

	imin=0; imax=34; imaxS=25;  imaxM=22; 	# should be imaxS=25.2 and imaxM=22.2 but it's to make the program faster
	di=1;
	Ni=Int(imax*(1/di))+1;
	NiS=Int(imaxS*(1/di))+1;
	NiM=Int(imaxM*(1/di))+1;
	i=collect(imin:di:imax);
	isympt=5; 					# should be 5.2
	Tsympt=Int(isympt*(1/di))+1;
	NiSympt=Int((imaxS-isympt)*(1/di))+1;

		# Parameters

	munat=zeros(Na);
	muadd=zeros(Na);
	gammaIdir=zeros(Na,Ni);
	gammaIindir=zeros(Na,Ni);

	betaS=zeros(1,Ni);
	betaM=zeros(1,Ni);
	betaA=zeros(1,Ni);
	nu=0.1;

	xiS=ones(Ni);					# case isolation of severe
	xiS[Tsympt:end]=exp.(-log(10)*(collect(isympt:di:imax)-isympt*ones(Ni-Tsympt+1)));
	xiM=ones(Ni);					# case isolation of mild
	xiM[Tsympt:end]=exp.(-log(2)*(collect(isympt:di:imax)-isympt*ones(Ni-Tsympt+1)));

	hS=ones(Na,Ni);   				# recovery rate for severe
	hM=ones(Na,Ni);					# recovery rate for mild (and paucisymptomatic)
	hS[:,1:NiS]=zeros(Na,NiS);
	hM[:,1:NiM]=zeros(Na,NiM);

	prophosp=zeros(Na);			# proportion of symptomatic

	if country=="France"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273,0.273,0.273]; # (data : Fergusson (2020))
		munat=[1-Nroot(365,1-0.00023),1-Nroot(365,1-0.00007),1-Nroot(365,1-0.00008),1-Nroot(365,1-0.00021),1-Nroot(365,1-0.00040), 1-Nroot(365,1-0.00045),1-Nroot(365,1-0.00057),1-Nroot(365,1-0.00081),1-Nroot(365,1-0.00125),1-Nroot(365,1-0.00204), 1-Nroot(365,1-0.00332), 1-Nroot(365,1-0.0052), 1-Nroot(365,1-0.0078), 1-Nroot(365,1-0.0106), 1-Nroot(365,1-0.0157), 1-Nroot(365,1-0.0217), 1-Nroot(365,1-0.0506),1-Nroot(365,1-0.0686),1-Nroot(365,1-0.151), 1-Nroot(365,1-0.231)];
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245),1-Nroot(NiSympt,1-0.3545),1-Nroot(NiSympt,1-0.3845)]; 		# I had to smooth a little bit the data in the end, because the last data was 80 years old and more

		if p==0.5
		 	alpha=1.378*10^(-7);				# with this apprroximation we need alpha=1.378*10^(-7) instead of 1.349*10^(-7) to have R0=3.3
		elseif p==0.2
			alpha=0.947*10^(-7);
		elseif p==0;
			alpha=0.783*10^(-7);
		else
			alpha=1.378*10^(-7);
			println("Because p is neither 0, 0.2 nor 0.5, the value of alpha must be changed to have R0=3.3.")
		end

	elseif country=="Burkina"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273]; # (data : Fergusson (2020))
		munat=[1-Nroot(365,1-0.0138),1-Nroot(365,1-0.0031),1-Nroot(365,1-0.0017),1-Nroot(365,1-0.0024),1-Nroot(365,1-0.0034), 1-Nroot(365,1-0.0041),1-Nroot(365,1-0.0048),1-Nroot(365,1-0.0059),1-Nroot(365,1-0.0071),1-Nroot(365,1-0.0092), 1-Nroot(365,1-0.0122), 1-Nroot(365,1-0.0174), 1-Nroot(365,1-0.0252), 1-Nroot(365,1-0.0385), 1-Nroot(365,1-0.0578), 1-Nroot(365,1-0.0847), 1-Nroot(365,1-0.1262),1-Nroot(365,1-0.2032)];
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245)];
		alpha=1.37*10^(-7);					# to have R0=3.3
	elseif country=="Vietnam"
		prophosp=[0.001,0.001,0.003,0.003,0.012,0.012,0.032,0.032,0.049,0.049,0.102,0.102,0.166,0.166,0.243,0.243,0.273,0.273,0.273,0.273]; # (data : Fergusson (2020))
		munat=(1-Nroot(365,1-6.317/1000))*ones(Na);
		gammaI_age=[1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.025),1-Nroot(NiSympt,1-0.0295),1-Nroot(NiSympt,1-0.0335),1-Nroot(NiSympt,1-0.0511),1-Nroot(NiSympt,1-0.0711),1-Nroot(NiSympt,1-0.117),1-Nroot(NiSympt,1-0.157),1-Nroot(NiSympt,1-0.196),1-Nroot(NiSympt,1-0.236),1-Nroot(NiSympt,1-0.2945),1-Nroot(NiSympt,1-0.3245),1-Nroot(NiSympt,1-0.3545),1-Nroot(NiSympt,1-0.3845)]; 		# I had to smooth a little bit the data in the end, because the last data was 80 years old and more
		alpha=0.616*10^(-7);				# to have R0=3.3 (Pasteur paper and reports of ETE team, Montpellier)
	end

	gammaI_infec=zeros(Ni); 				# age of infection distribution of deaths due to COVID-19
	gammaI_infec[Tsympt:NiS]=ones(NiSympt);	# death due to COVID-19 only occur between i_sympt and i_max for severe infections
	gammaIdir=gammaI_age*gammaI_infec';

	betak=3;betalambda=5/(log(2)^(1/betak));

	betaAA=(alpha)*(betak/betalambda)*((i/betalambda).^(betak-1)).*exp.(-(i/betalambda).^betak); # transmission probabilities
	betaS=betaAA.*xiS;
    betaM=betaAA.*xiM;
	betaA=betaAA*nu;

		# Social contact matrix for France

	betacontact_age=ones(Na,Na);

	if country=="France"
		betacontact_age=[3.805	 1.106	 0.412	 0.255	 0.342	 0.753	 1.209	 1.042	 0.496	 0.31	 0.307	 0.272	 0.186	 0.128	 0.088	 0.054	 0.054	 0.054	 0.054	 0.054;
			1.062	 5.033	 1.011	 0.273	 0.164	 0.415	 0.921	 1.138	 0.845	 0.331	 0.207	 0.185	 0.177	 0.112	 0.058	 0.053	 0.053	 0.053	 0.053	 0.053;
			0.238	 1.538	 6.986	 0.846	 0.287	 0.244	 0.403	 0.827	 1.081	 0.57	 0.266	 0.147	 0.102	 0.088	 0.07	 0.07	 0.07	 0.07	 0.07	 0.07;
			0.124	 0.308	 2.301	 7.832	 1.36	 0.651	 0.531	 0.762	 1.025	 1.049	 0.479	 0.164	 0.072	 0.053	 0.035	 0.026	 0.026	 0.026	 0.026	 0.026;
			0.2	 0.168	 0.256	 2.165	 3.934	 1.708	 1.168	 1.133	 0.97	 1.27	 0.813	 0.327	 0.084	 0.046	 0.054	 0.052	 0.052	 0.052	 0.052	 0.052;
			0.533	 0.249	 0.144	 0.731	 1.944	 3.458	 1.777	 1.471	 1.263	 0.997	 0.961	 0.357	 0.113	 0.052	 0.031	 0.024	 0.024	 0.024	 0.024	 0.024;
			0.722	 0.848	 0.53	 0.423	 1.005	 1.72	 2.919	 1.894	 1.457	 1.135	 0.81	 0.405	 0.184	 0.085	 0.047	 0.053	 0.053	 0.053	 0.053	 0.053;
			0.704	 1.072	 0.837	 0.683	 0.765	 1.407	 1.76	 3.211	 2.098	 1.344	 0.92	 0.364	 0.26	 0.152	 0.098	 0.046	 0.046	 0.046	 0.046	 0.046;
			0.31	 0.657	 0.983	 1.109	 0.924	 1.252	 1.591	 1.877	 2.954	 1.682	 1.113	 0.273	 0.184	 0.111	 0.087	 0.048	 0.048	 0.048	 0.048	 0.048;
			0.406	 0.467	 0.612	 1.537	 0.945	 0.994	 1.259	 1.465	 1.573	 2.196	 1.086	 0.341	 0.151	 0.086	 0.083	 0.082	 0.082	 0.082	 0.082	 0.082;
			0.252	 0.595	 0.829	 1.202	 1.019	 1.291	 1.179	 1.215	 1.594	 1.743	 1.918	 0.649	 0.267	 0.119	 0.088	 0.09	 0.09	 0.09	 0.09	 0.09;
			0.542	 0.671	 0.574	 0.749	 0.614	 0.932	 0.918	 0.748	 0.793	 0.681	 0.923	 1.462	 0.513	 0.218	 0.12	 0.095	 0.095	 0.095	 0.095	 0.095;
			0.396	 0.363	 0.257	 0.361	 0.339	 0.508	 0.639	 0.694	 0.523	 0.416	 0.441	 0.675	 1.406	 0.426	 0.263	 0.129	 0.129	 0.129	 0.129	 0.129;
			0.198	 0.313	 0.257	 0.157	 0.212	 0.307	 0.45	 0.454	 0.425	 0.302	 0.324	 0.44	 0.518	 1.003	 0.276	 0.139	 0.139	 0.139	 0.139	 0.139;
			0.105	 0.294	 0.308	 0.355	 0.156	 0.245	 0.253	 0.477	 0.515	 0.434	 0.339	 0.338	 0.66	 0.619	 0.974	 0.326	 0.326	 0.326	 0.326	 0.326;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643;
			0.244	 0.317	 0.452	 0.36	 0.154	 0.188	 0.305	 0.382	 0.451	 0.443	 0.539	 0.344	 0.28	 0.367	 0.368	 0.643	 0.643	 0.643	 0.643	 0.643];
	elseif country=="Burkina"
		betacontact_age[1:3,:]=[6.415 3.639 1.748 0.919 1.198 1.574 1.469 1.158 0.661 0.398 0.433 0.377 0.251 0.163 0.079 0.036 0.036 0.036;3.477 18.279 4.124 1.077 0.652 1.292 1.3506 1.237 0.911 0.452 0.356 0.345 0.283 0.160 0.068 0.037 0.037 0.037;1.177 5.907 17.935 1.868 0.832 0.717 0.830 0.862 0.851 0.500 0.314 0.169 0.1295 0.113 0.068 0.033 0.033 0.033];
		betacontact_age[4:6,:]=[0.634 1.532 5.808 13.816 2.361 1.144 0.784 0.923 0.852 0.703 0.440 0.228 0.141 0.064 0.031 0.015 0.015 0.015;0.849 0.858 0.911 4.087 4.586 2.153 1.273 1.049 0.824 0.768 0.581 0.346 0.179 0.038 0.029 0.017 0.017 0.017;1.172 0.782 0.525 1.225 2.435 2.846 1.660 1.282 1.038 0.797 0.725 0.407 0.221 0.038 0.013 0.010 0.010 0.010];
		betacontact_age[7:9,:]=[1.125 1.647 1.311 0.664 1.163 1.753 1.842 1.502 1.175 0.902 0.657 0.471 0.248 0.043 0.021 0.013 0.013 0.013;0.900 1.390 1.052 0.789 0.725 1.249 1.356 1.675 1.461 0.962 0.744 0.385 0.200 0.063 0.028 0.009 0.009 0.0091;0.738 1.218 1.200 1.010 0.833 1.094 1.302 1.365 1.516 1.168 0.891 0.366 0.255 0.077 0.029 0.013 0.013 0.013];
		betacontact_age[10:12,:]=[0.464 0.922 0.951 1.188 0.690 0.873 1.009 1.117 1.070 1.035 0.781 0.448 0.186 0.044 0.027 0.022 0.022 0.022;0.565 0.920 1.185 1.005 0.724 1.029 1.013 0.980 1.238 1.241 1.001 0.639 0.234 0.040 0.030 0.022 0.022 0.022;1.290 1.653 1.259 0.956 0.700 1.103 1.183 0.908 0.985 0.923 0.928 0.627 0.346 0.108 0.033 0.038 0.038 0.038];
		betacontact_age[13:15,:]=[0.827 0.941 0.750 0.525 0.443 0.629 0.684 0.775 0.693 0.622 0.575 0.506 0.219 0.098 0.039 0.010 0.010 0.010;0.575 0.935 0.720 0.257 0.232 0.246 0.335 0.357 0.284 0.163 0.189 0.174 0.126 0.083 0.059 0.017 0.017 0.017;0.244 0.792 0.6343 0.393 0.103 0.205 0.141 0.258 0.241 0.238 0.243 0.144 0.154 0.099 0.061 0.061 0.061 0.061];
		betacontact_age[16:18,:]=[0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022;0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022;0.311 0.453 0.559 0.369 0.095 0.106 0.124 0.224 0.184 0.199 0.217 0.149 0.061 0.090 0.045 0.022 0.022 0.022];
	elseif country=="Vietnam"
		betacontact_age[1:3,:]=[3.424 1.396 0.631 0.441 0.655 1.08 1.26 1.078 0.624 0.368 0.39 0.303 0.168 0.103 0.062 0.037 0.037 0.037 0.037 0.037;1.316 6.409 1.498 0.479 0.306 0.767 1.069 1.121 0.92 0.401 0.27 0.232 0.16 0.093 0.046 0.037 0.037 0.037 0.037 0.037;0.406 2.212 8.116 1.164 0.506 0.465 0.678 0.926 1.088 0.606 0.335 0.156 0.086 0.075 0.058 0.046 0.046 0.046 0.046 0.046];
		betacontact_age[4:6,:]=[0.258 0.608 3.35 11.788 2.213 1.134 0.848 1.088 1.242 1.139 0.614 0.253 0.1 0.051 0.03 0.022 0.022 0.022 0.022 0.022;0.414 0.345 0.485 3.478 5.118 2.538 1.595 1.409 1.168 1.284 0.908 0.471 0.147 0.038 0.041 0.035 0.035 0.035 0.035 0.035;0.856 0.492 0.329 1.261 2.945 4.328 2.415 1.858 1.52 1.23 1.165 0.616 0.223 0.054 0.024 0.022 0.022 0.022 0.022 0.022];
		betacontact_age[7:9,:]=[0.894 1.215 0.953 0.684 1.419 2.386 3.173 2.303 1.72 1.307 0.981 0.654 0.269 0.072 0.039 0.033 0.033 0.033 0.033 0.033;0.857 1.315 1.14 1.023 0.994 1.809 2.159 3.085 2.335 1.483 1.054 0.537 0.256 0.118 0.066 0.026 0.026 0.026 0.026 0.026;0.523 0.952 1.216 1.406 1.192 1.544 1.922 2.135 2.726 1.805 1.253 0.455 0.221 0.094 0.058 0.031 0.031 0.031 0.031 0.031];
		betacontact_age[10:12,:]=[0.462 0.693 0.856 1.766 1.116 1.256 1.459 1.652 1.697 1.892 1.177 0.543 0.164 0.056 0.052 0.053 0.053 0.053 0.053 0.053;0.467 0.767 1.072 1.487 1.282 1.659 1.466 1.394 1.8 1.864 1.663 0.879 0.24 0.064 0.052 0.055 0.055 0.055 0.055 0.055;0.772 0.894 0.783 0.992 0.888 1.42 1.464 1.129 1.225 1.026 1.156 0.992 0.356 0.114 0.05 0.044 0.044 0.044 0.044 0.044];
		betacontact_age[13:15,:]=[0.571 0.533 0.394 0.481 0.452 0.69 0.786 0.818 0.633 0.527 0.484 0.5 0.392 0.167 0.077 0.03 0.03 0.03 0.03 0.03;0.304 0.46 0.372 0.204 0.227 0.301 0.447 0.463 0.374 0.19 0.19 0.218 0.193 0.253 0.093 0.034 0.034 0.034 0.034 0.034;0.155 0.444 0.389 0.364 0.136 0.256 0.238 0.429 0.421 0.326 0.247 0.173 0.223 0.183 0.237 0.111 0.111 0.111 0.111 0.111];
		betacontact_age[16:end,:]=[0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131;0.293 0.409 0.59 0.505 0.167 0.179 0.25 0.398 0.408 0.384 0.361 0.193 0.098 0.147 0.117 0.131 0.131 0.131 0.131 0.131];
	end

		# Initial conditions

	S=zeros(nbr,Na);		# Susceptible
	Is=zeros(nbr,Na,Ni);	# Severe infected
	Im=zeros(nbr,Na,Ni);	# Mild infected
	Ip=zeros(nbr,Na,Ni);	# Paucisymptomatic infected

	H=zeros(nbr);			# hospitalised
	R=zeros(nbr,Na);		# recovered
	Mindir=zeros(nbr,Na);	# deaths directly due to COVID-19
	Mdir=zeros(nbr,Na);		# deaths indirectly due to COVID-19
	Ms=zeros(nbr,Na);		# natural mortality

	if country=="France"
			prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.03,0.02,0.01]; # age distribution of hospitalised (mean data between reports of Santé Publique France and CDC
			Is[1,:,1:NiS]=(130/(NiS*di))*prop_age*ones(NiS)'; 	# 130 is the number of severe infected in France at March, 1st, uniformly distributed in age since infection then distributed in age
			Im[1,:,1:NiS]=(130/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
			Ip[1,:,1:NiS]=(130/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
			S[1,:]=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
		elseif country=="Burkina"
			prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.06]; #
			Is[1,:,1:NiS]=(288/(NiS*di))*prop_age*ones(NiS)'; 	# 288 is the number of cases in Burkina, on April, 1st (WHO ref)
			Im[1,:,1:NiS]=(288/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
			Ip[1,:,1:NiS]=(288/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
			S[1,:]=18450494/100*[13.2,12.6,11.6,12.9,11.5,9.3,7.2,5.4,4.3,3.1,2.5,1.8,1.4,0.9,0.7,0.4,0.3,0.2]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
		elseif country=="Vietnam"
			prop_age=[0.005,0.005,0.01,0.03,0.05,0.05,0.05,0.06,0.06,0.09,0.09,0.08,0.08,0.07,0.07,0.07,0.07,0.03,0.02,0.01]; # age distribution of hospitalised (mean data between reports of Santé Publique France and CDC
			Is[1,:,1:NiS]=(217/(NiS*di))*prop_age*ones(NiS)'; 	# 130 is the number of severe infected in France at March, 1st, uniformly distributed in age since infection then distributed in age
			Im[1,:,1:NiS]=(217/(NiS*di))*(prop_age./prophosp).*(ones(Na)-prophosp)*ones(NiS)'; 	# we deduce the mild infected
			Ip[1,:,1:NiS]=(217/(NiS*di))*(prop_age./prophosp)*(p/(1-p))*ones(NiS)';				# and the paucisymptomatic
			S[1,:]=[4539031,7525542,6976064,6474991,7145151,8741274,8345522,7631483,6941367,6442420,5774836,5136526,4149565,2772903,1526235,1113819,873832,606471,268159,124992]'-di*sum(Is[1,:,:],dims=2)[:,1]'-di*sum(Im[1,:,:],dims=2)[:,1]'-di*sum(Ip[1,:,:],dims=2)[:,1]';	# population age distributed at the beginning of the epidemic
		end

	H[1]=di*sum(Is[1,:,Tsympt:end]);

	muadd=(0.01*munat)/(1+99*exp(-10*(((H[1])/(Hsat))-1)));
	gammaIindir=(gammaIdir)/(1+99*exp(-10*(((H[1])/(Hsat))-1)));

	Mdir[1,:]=sum(gammaIdir.*Is[1,:,:],dims=2)[:,1];
	Mindir[1,:]=sum(gammaIindir.*Is[1,:,:],dims=2)[:,1]+muadd.*S[1,:]+di*sum(muadd.*Is[1,:,:],dims=2)[:,1]+di*sum(muadd.*Im[1,:,:],dims=2)[:,1]+di*sum(muadd.*Ip[1,:,:],dims=2)[:,1]+muadd.*R[1,:];
	Ms[1,:]=muadd.*S[1,:]+munat.*S[1,:];

			# R0 computation

 	R0=0;
	# piS=zeros(Na,Ni);
	# piM=zeros(Na,Ni);
	# pIp=zeros(Na,Ni);
	# omega=zeros(Na,Ni);
	#
	# for aa=1:Na
	# 	piS[aa,:]=exp.(-munat[aa]*i-(i-imaxS*ones(Ni)).*Charac.(i,imaxS,imax+di)-gammaIdir[aa,:].*(i-isympt*ones(Ni)).*Charac.(i,isympt,imaxS));
	# 	piM[aa,:]=exp.(-munat[aa]*i-(i-imaxM*ones(Ni)).*Charac.(i,imaxM,imax+di));
	# 	pIp[aa,:]=exp.(-munat[aa]*i-(i-imaxM*ones(Ni)).*Charac.(i,imaxM,imax+di));
	# 	omega[aa,:]=(1-p)*prophosp[aa]*betaS.*piS[aa,:]+(1-p)*(1-prophosp[aa])*betaM.*piM[aa,:]+p*betaA.*pIp[aa,:];
	# end
	#
	# Omega=zeros(Na)
	# Omega=di*sum(omega,dims=2)[:,1];
	# Omegamat=hcat(fill.(Omega,Na)...)';
	#
	# S0mat=hcat(fill.(S[1,:],Na)...)';
	#
	# R0mat=betacontact_age.*Omegamat.*S0mat;
	# R0=maximum(abs.(eigvals(R0mat)));		# need package LinearAlgebra to compute eigenvalues

		# Confinement

	betacontact=ones(Na,Na);

		### Algorithm ###

	S=convert(Array{Float64,2},S);
	Is=convert(Array{Float64,3},Is);
	Im=convert(Array{Float64,3},Im);
	Ip=convert(Array{Float64,3},Ip);

	Lambda=zeros(nbr,Na);			# force of infection (age dependant)
	Lambda2=zeros(nbr,Na);			# force of infection without control measures (to compute optimal later)

	for tt=1:(nbr-1)

		betacontact=betacontact_age.*(ones(Na,Na)-hcat(fill.(c[tt,:],Na)...)');

			# We compute mortality rates, depending on the number of hospitalisations

		muadd=(0.01*munat)/(1+99*exp.(-10*(((H[tt])/(Hsat))-1)));
		gammaIindir=(gammaIdir)/(1+99*exp(-10*(((H[tt])/(Hsat))-1)));

			# We continue the algorithm ##

		Lambda1=di*Is[tt,:,:]*betaS+di*Im[tt,:,:]*betaM+di*Ip[tt,:,:]*betaA; 	# we compute the force of infection according to the age
		Lambda[tt,:]=betacontact*Lambda1;
		Lambda2[tt,:]=betacontact_age*Lambda1;

		S[tt+1,:]=(S[tt,:])./(ones(Na)+dt*munat+dt*muadd+dt*Lambda[tt,:]); 		# we compute new density of susceptible
		Is[tt+1,:,1]=(Is[tt,:,1]+(1-p)*(dt/di)*prophosp.*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hS[:,1]+(dt/di)*gammaIdir[:,1]+(dt/di)*gammaIindir[:,1]);
		Im[tt+1,:,1]=(Im[tt,:,1]+(1-p)*(dt/di)*(ones(Na)-prophosp).*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,1]);
		Ip[tt+1,:,1]=(Im[tt,:,1]+p*(dt/di)*Lambda[tt,:].*S[tt+1,:])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,1]);

			# Loop #

		for ii=1:(Ni-1)
			Is[tt+1,:,ii+1]=(Is[tt,:,ii+1]+(dt/di)*Is[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hS[:,ii+1]+(dt/di)*gammaIdir[:,ii+1]+(dt/di)*gammaIindir[:,ii+1]);
			Im[tt+1,:,ii+1]=(Im[tt,:,ii+1]+(dt/di)*Im[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,ii+1]);
			Ip[tt+1,:,ii+1]=(Ip[tt,:,ii+1]+(dt/di)*Ip[tt+1,:,ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+dt*hM[:,ii+1]);
		end

		H[tt+1]=di*sum(Is[tt,:,Tsympt:end]);

		Mdir[tt+1,:]=Mdir[tt,:]+dt*sum(gammaIdir.*Is[tt+1,:,:],dims=2)[:,1];
		Mindir[tt+1,:]=Mindir[tt,:]+dt*sum(gammaIindir.*Is[tt+1,:,:],dims=2)[:,1]+dt*muadd.*S[tt+1,:]+dt*di*sum(muadd.*Is[tt+1,:,:],dims=2)[:,1]+dt*di*sum(muadd.*Im[tt+1,:,:],dims=2)[:,1]+dt*di*sum(muadd.*Ip[tt+1,:,:],dims=2)[:,1]+dt*muadd.*R[tt+1,:];
		R[tt+1,:]=(R[tt,:]+dt*Is[tt,:,Ni]+dt*Im[tt,:,Ni]+dt*Ip[tt,:,Ni]+dt*di*sum(hS.*Is[tt+1,:,:],dims=2)[:,1]+dt*di*sum(hM.*Im[tt+1,:,:],dims=2)[:,1]+dt*di*sum(hM.*Ip[tt+1,:,:],dims=2)[:,1])./(ones(Na)+dt*munat+dt*muadd);
		Ms[tt,:]=dt*muadd.*S[tt,:]+dt*munat.*S[tt,:];

	end

	###### Adjoint system ######

	c_hat=zeros(nbr,Na);

	zS=zeros(nbr,Na);
	zIs=zeros(nbr,Na,Ni);
	zIm=zeros(nbr,Na,Ni);
	zIp=zeros(nbr,Na,Ni);
	zR=zeros(nbr,Na);
	zeta1=zeros(nbr,Na);
	zeta2=zeros(nbr,Na);
	zeta3=zeros(nbr,Na);

	for tt=1:(nbr-1)

		muadd=(0.01*munat)/(1+99*exp.(-10*(((H[nbr-tt+1])/(Hsat))-1)));
		gammaIindir=(gammaIdir)/(1+99*exp(-10*(((H[nbr-tt+1])/(Hsat))-1)));

		muadd_der=(0.01*munat*10*(99/Hsat)*exp.(-10*(((H[nbr-tt+1])/(Hsat))-1)))/((1+99*exp(-10*(((H[nbr-tt+1])/(Hsat))-1)))^2);
		gammaIindir_der=(gammaIdir*10*(99/Hsat)*exp.(-10*(((H[nbr-tt+1])/(Hsat))-1)))/((1+99*exp(-10*(((H[nbr-tt+1])/(Hsat))-1)))^2);
		zeta1[nbr-tt,:]=di*muadd_der.*(sum(Im[nbr-tt,:,:].*(ones(Na,Ni)-zIm[nbr-tt+1,:,:]),dims=2)[:,1]+sum(Ip[nbr-tt,:,:].*(ones(Na,Ni)-zIp[nbr-tt+1,:,:]),dims=2)[:,1]+sum(Is[nbr-tt,:,:].*(ones(Na,Ni)-zIs[nbr-tt+1,:,:]),dims=2)[:,1]+S[nbr-tt,:].*(ones(Na)-zS[nbr-tt+1,:])+R[nbr-tt,:].*(ones(Na)-zR[nbr-tt+1,:]))+sum(gammaIindir_der.*Is[nbr-tt,:,:].*(ones(Na,Ni)-zIs[nbr-tt+1,:,:]),dims=2)[:,1];
		zeta2[nbr-tt,:]=(ones(Na)-c[nbr-tt,:]).*(((1-p)*(ones(Na)-prophosp).*zIm[nbr-tt+1,:,1]+(1-p)*prophosp.*zIs[nbr-tt+1,:,1]+p*zIp[nbr-tt+1,:,1])-zS[nbr-tt+1,:]);
		zeta3[nbr-tt,:]=zR[nbr-tt,:];
		zS[nbr-tt,:]=(zS[nbr-tt+1,:]+dt*muadd+di*dt*betacontact_age*(Is[nbr-tt+1,:,:]*betaS+Im[nbr-tt+1,:,:]*betaM+Ip[nbr-tt+1,:,:]*betaA).*zeta2[nbr-tt,:])./(ones(Na)+dt*munat+dt*muadd);
		zR[nbr-tt,:]=(zR[nbr-tt+1,:]+dt*muadd)./(ones(Na)+dt*munat+dt*muadd);
		zIs[nbr-tt,:,Ni]=(zIs[nbr-tt+1,:,Ni]+dt*muadd+dt*hS[:,Ni]+(dt/di)*gammaIdir[:,Ni]+(dt/di)*gammaIindir[:,Ni]+dt*Charac((Ni-1)*di,isympt,imaxS)*zeta1[nbr-tt,:]+dt*(betaS[Ni]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hS[:,Ni])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+(dt/di)*gammaIdir[:,Ni]+(dt/di)*gammaIindir[:,Ni]);
		zIm[nbr-tt,:,Ni]=(zIm[nbr-tt+1,:,Ni]+dt*muadd+dt*hM[:,Ni]+dt*(betaM[Ni]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hM[:,Ni])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd);
		zIp[nbr-tt,:,Ni]=(zIp[nbr-tt+1,:,Ni]+dt*muadd+dt*hM[:,Ni]+dt*(betaA[Ni]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hM[:,Ni])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd);

		for ii=1:(Ni-1)
			zIs[nbr-tt,:,Ni-ii]=(zIs[nbr-tt+1,:,Ni-ii]+(dt/di)*zIs[nbr-tt,:,Ni-ii+1]+dt*muadd+dt*hS[:,Ni-ii]+(dt/di)*gammaIdir[:,Ni-ii]+(dt/di)*gammaIindir[:,Ni-ii]+dt*Charac((Ni-ii-1)*di,isympt,imaxS)*zeta1[nbr-tt,:]+dt*(betaS[Ni-ii]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hS[:,Ni-ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd+(dt/di)*gammaIdir[:,Ni-ii]+(dt/di)*gammaIindir[:,Ni-ii]);
			zIm[nbr-tt,:,Ni-ii]=(zIm[nbr-tt+1,:,Ni-ii]+(dt/di)*zIm[nbr-tt,:,Ni-ii+1]+dt*muadd+dt*hM[:,Ni-ii]+dt*(betaM[Ni-ii]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hM[:,Ni-ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd);
			zIp[nbr-tt,:,Ni-ii]=(zIp[nbr-tt+1,:,Ni-ii]+(dt/di)*zIp[nbr-tt,:,Ni-ii+1]+dt*muadd+dt*hM[:,Ni-ii]+dt*(betaA[Ni-ii]*betacontact_age*zeta2[nbr-tt,:]).*S[nbr-tt,:]+dt*zeta3[nbr-tt,:].*hM[:,Ni-ii])./(ones(Na)+(dt/di)*ones(Na)+dt*munat+dt*muadd);
		end

		c_hat[nbr-tt,:]=(S[nbr-tt,:].*Lambda2[nbr-tt,:].*((1-p)*(ones(Na)-prophosp).*zIm[nbr-tt,:,1]+(1-p)*prophosp.*zIs[nbr-tt,:,1]+p*zIp[nbr-tt,:,1])-zS[nbr-tt,:])./(2*B);

	end

	prix=sum(Mdir[end,:])+sum(Mindir[end,:]);

	c_temp=sum(c_hat[:,1])*dt;


	#if cseuil==0
		c_hat=max.(min.(c_hat,0.95),0);
		#c_hat[:,1:6]=zeros(nbr,6);
		#c_hat[:,1:10]=zeros(nbr,10);
	#else
	#	c_sum=dt*sum(((c_hat).^2)*B);
	#	c_hat=c_hat*cseuil/c_sum;
	#	c_hat=max.(min.(c_hat,1),0);
	#end

	return c_hat,prix,S,Is,Im,Ip,R,Lambda,Lambda2,zS,zIs,zIm,zIp,zR,Mdir,Mindir,R0;
end

		##### Optimal algorithm #####

function CovidOptim(N,Bstar,ctest,p,country) 	## N is the number of iterations

 ## [400,430,460],Bstar,0

	T=365;
	dt=0.2;	t=collect(0:dt:T);nbr=length(t);
	amin=0;

	if country=="France"
		amax=19;
		S=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655];
	elseif country=="Burkina"
		amax=17;
		S=18450494/100*[13.2,12.6,11.6,12.9,11.5,9.3,7.2,5.4,4.3,3.1,2.5,1.8,1.4,0.9,0.7,0.4,0.3,0.2];
	elseif country=="Vietnam"
		amax=19;
		S=[4539031,7525542,6976064,6474991,7145151,8741274,8345522,7631483,6941367,6442420,5774836,5136526,4149565,2772903,1526235,1113819,873832,606471,268159,124992];
	end

	a=collect(0:1:amax);Na=length(a);

	# We put uniform cost for control measures (uniform by individuals and not by class of age)

	sumS=sum(S);
	B=Bstar*S/sumS;

	size_N=length(N);
	Nmax=maximum(N);
	Nmax_err0=Nmax;

	err=zeros(Nmax); 		# error between two iterations

	if length(ctest)==1		# so that if we put a uniform control (c=cst), then c becomes a constant matrix
		ctest=ctest*ones(nbr,Na);
	end

	c=zeros(nbr,Na,Nmax+1);
	prix1=zeros(Nmax+1);	# price to pay (deaths + cost of measures control)
	prix2=zeros(Nmax+1);
	c[:,:,1]=ctest;			# initialisation

	## We start the first iterations without convex combinations

	cvx=0;
	ii=1;
	convexe=cvx*ones(nbr,Na);

	while ii<3

		cc,prix1[ii]=Covidadj(T,c[:,:,ii],B,p,country);
		prix2[ii]=dt*sum((c[:,:,ii].^2)*B);

		c[:,:,ii+1]=convexe.*cc+(ones(nbr,Na)-convexe).*c[:,:,ii];
		err[ii]=sum(abs.(c[:,:,ii+1]-c[:,:,ii]))/sum(abs.(c[:,:,ii+1]));

		ii=ii+1;

	end

	## We take the best iteration among the three firsts, and we continue with convex iterations

	f=argmin((prix1+prix2)[1:2]);
	c[:,:,3]=c[:,:,f];

	convexe=hcat(0.2*fill.((collect(0:dt:T)./(T)),Na)...)';

	## It's a convex combination in time, but also in age.

	while ii<N[1]+1

		vectt=zeros(Na);

		for ll=1:Na
			kk=1;
			vectt[Na-ll+1]=1;
			convexe=(collect(0:dt:T)/T)*(vectt');

			while (kk<21)&(ii<N[1]+1)

				cc,prix1[ii]=Covidadj(T,c[:,:,ii],B,p,country);
				prix2[ii]=dt*sum((c[:,:,ii].^2)*B);
				c[:,:,ii+1]=convexe.*cc+(ones(nbr,Na)-convexe).*c[:,:,ii];
				err[ii]=sum(abs.(c[:,:,ii+1]-c[:,:,ii]))/sum(abs.(c[:,:,ii+1]));

				if 0==ii%20
					println("ii=$(ii)")	# to know where we are in the algorithm process
				end

				ii=ii+1;
				kk=kk+1;
			end
			ff=argmin((prix1+prix2)[1:(ii-1)]);
			c[:,:,ii]=c[:,:,ff];
		end
	end

		## Now we take the mean of the 10 last iterations, to see if we get better control

	if (size_N>1)

		if N[2]>N[1]

			while (ii<N[2]+1)&(err[ii-1]!=0)

				cc,prix1[ii]=Covidadj(T,c[:,:,ii],B,p,country);
				prix2[ii]=dt*sum((c[:,:,ii].^2)*B);
				c[:,:,ii+1]=(sum(c[:,:,ii-8:ii],dims=3)[:,:,1]+cc)/10;
				err[ii]=sum(abs.(c[:,:,ii+1]-c[:,:,ii]))/sum(abs.(c[:,:,ii+1]));

				if 0==ii%20
					println("ii=$(ii)")
				end

				if err[ii-1]==0
					Nmax_err0=ii;
				end

				ii=ii+1;

			end
		end
	end

		## Now we continue the convex combinations (convex only in time)

	convexe=hcat(0.2*fill.((collect(0:dt:T)./(T)),Na)...)';

	if (size_N>2)
		if N[3]>N[2]

			ff=argmin((prix1+prix2)[1:(ii-1)]);
			c[:,:,ii]=c[:,:,ff];

			while (ii<N[3]+1)&(err[ii-1]!=0)

				cc,prix1[ii]=Covidadj(T,c[:,:,ii],B,p,country);
				prix2[ii]=dt*sum((c[:,:,ii].^2)*B);
				c[:,:,ii+1]=convexe.*cc+(ones(nbr,Na)-convexe).*c[:,:,ii];
				err[ii]=sum(abs.(c[:,:,ii+1]-c[:,:,ii]))/sum(abs.(c[:,:,ii+1]));

				if 0==ii%20
					println("ii=$(ii)")
				end

				ii=ii+1;
			end
		end
	end

	mini=argmin((prix1+prix2)[1:(Nmax_err0-1)]);
	pp3D=plot(t,a*5,c[:,:,mini]',st=:surface,xlabel="Time (days)",ylabel="Age (years)",title="Optimal control for \$B^{\\ast}\$=$(Bstar)",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,clims=(0,1),xticks=collect(0:50:T),yticks=collect(0:10:100),zticks=(collect(0:0.2:1)),zlims=(0,1),ztickfont=font(30),zguidefontsize=34,colorbar=false);
	ppciel=plot(t,a*5,c[:,:,mini]',st=:contourf,xlabel="Time (days)",ylabel="Age (years)",title="Optimal control for \$B^{\\ast}\$=$(Bstar)",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,clims=(0,1),xticks=collect(0:50:T),yticks=collect(0:20:100));
	pp=plot(pp3D,ppciel,layout=@layout([pp3D{0.6w} ppciel]));

	return c,err,ppciel,prix1[1:(end-1)],prix2[1:(end-1)];	# return the optimal control c, the error err and some figures.

end

		##### Optimal figures - Figures 6, S2, S3 when p=0.5 and for each country - and Figure S5 when p=0 or p=0.2, for France and B=1000 #####

function CovidOptimFig(N,Bstar,p,country)

	### Values

	### N=[400,430,460];
	### Bstar=[10000,1000,100];

	### For p=0.5

	### argmin=[457,225,432] respectively, min=[14808,92100,325888];
	### prix1=[230,42934,143662] (605939 when c=0) and prix2=[14578,49166,182224]

	if country=="France"
		Hsat=5000;
	elseif country=="Burkina"
		Hsat=11;
	elseif country=="Vietnam"
		Hsat=5932;
	end

	timecomput=time();	# Initial time
	size_Bstar=length(Bstar);
	T=365;
	ptot=zeros(size_Bstar+1);
	deathstot=zeros(size_Bstar+1);

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,0,Hsat,p,5.2,0.1,country);
	dt=0.2;	nbr=Int(T/dt+1);tt=ones(nbr);
	Na=length(a); aa=ones(Na);
	hosp_max=maximum(Hosp);
	couleurs=[:darkorange,:limegreen,:orchid];

	ptot[1]=prop_totale[end];
	deathstot[1]=sum(Mdir[end,:])+sum(Mindir[end,:]);

	p3=plot(t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),xlabel="Time (days)",ylabel="\$\\log_{10}\$(Number of cases)",title="Number of hospitalisations",titlefontsize=38,xtickfont=font(30),ytickfont=font(30),xguidefontsize=32,yguidefontsize=30,xticks=collect(0:30:365),color="blue",label="Without control",legendfont = font("Times new roman", 26),width=2);
	p4=plot(a*5,log10.(aa+(Mindir)[nbr,:]+(Mdir)[nbr,:]),xlabel="Age (years)",ylabel="\$\\log_{10}\$(Deaths)",title="Cumulative number of deaths at t=$(T)",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xticks=collect(0:20:100),color="blue",label="Without control",legendfont = font("Times new roman", 24),legend=:topleft,width=2);

	c_opt=zeros(nbr,Na,size_Bstar);
	pp_opt=[];

	for ii=1:size_Bstar

	c,err,ppciel,prix1,prix2=CovidOptim(N,Bstar[ii],0,p,country);
	pp_opt=push!(pp_opt,ppciel);
	it_min=argmin(prix1+prix2);
	c_opt[:,:,ii]=c[:,:,it_min];

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,c_opt[:,:,ii],Hsat,p,5.2,0.1,country);
	ptot[ii+1]=prop_totale[nbr];
	deathstot[ii+1]=sum(Mdir[end,:])+sum(Mindir[end,:]);

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="With \$B^{\\ast}\$=\$10^{$(5-ii)}\$",color=couleurs[ii],width=2);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="With \$B^{\\ast}\$=\$10^{$(5-ii)}\$",color=couleurs[ii],width=2);

	end

	p3=plot!(p3,t,log10.(tt+Hsat*ones(size(t))),linestyle=:dot,label="Healthcare capacity",color="darkorange",legend=:topright,width=2);

	println("The file took $(time()-timecomput) seconds to execute.")

	return pp_opt,p3,p4,c_opt,ptot,deathstot;

end

		##### Figures for different strategies - Figures 7, 8, 9 for each country - #####

function CovidOptimFigStrategies(T,country)

	timecomput=time();	# initial time

	Bstar=1000;
	TotalMort=zeros(4);
	ptot=zeros(4);

	if country=="France"
		Hsat=5000;
		amax=20;
	elseif country=="Burkina"
		Hsat=11;
		amax=18;
	elseif country=="Vietnam"
		Hsat=5932;
		amax=20;
	else
		println("You need to write correctly either France, Burkina or Vietnam")
		return
	end

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,0,Hsat,0.5,5.2,0.1,country);
	dt=0.2;	nbr=Int(T/dt+1);tt=ones(nbr);
	Na=length(a); aa=ones(Na);
	hosp_max=maximum(Hosp);
	couleurs=[:darkorange,:limegreen,:orchid];

	p3=plot(t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),xlabel="Time (days)",ylabel="\$\\log_{10}\$(Number of cases)",title="Number of hospitalisations",titlefontsize=38,xtickfont=font(30),ytickfont=font(30),xguidefontsize=32,yguidefontsize=30,xticks=collect(0:30:365),color="blue",label="Without control",legendfont = font("Times new roman", 26),width=2);

	p4=plot(a*5,log10.(aa+(Mindir)[nbr,:]+(Mdir)[nbr,:]),xlabel="Age (years)",ylabel="\$\\log_{10}\$(Deaths)",title="Cumulative number of deaths at t=$(T)",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xticks=collect(0:20:100),color="blue",label="Without control",legendfont = font("Times new roman", 28),width=2);

	TotalMort[1]=sum(Mdir[end,:])+sum(Mindir[end,:]);
	ptot[1]=prop_totale[end];

	# Proportions of infected for each strategy

	p5=plot(collect(0:5:((amax)*5)),[prop_age[end,:];prop_age[end,end]],xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xlabel="Age (years)",ylabel="Proportion",title="Proportion of infected at t=$T",color="blue",label="Without control",titlefontsize=38,ylims=(0,1),width=2,linetype=:steppost);
	p5=plot!(p5,xticks = (collect(0:5:T)),yticks=collect(0:0.1:1),legendfont = font("Times new roman", 28));
	p5=plot!(p5,left_margin=-1.9cm, right_margin=-2.5cm,bottom_margin=-0.9cm,top_margin=-0.3cm);
	p5=plot!(p5,collect(0:5:(amax)*5),prop_totale[end]*ones(Na),label="",linestyle=:dot,color="blue",width=2);

	# Optimal strategy

	c,err,ppciel,prix1,prix2=CovidOptim([400,430,460],Bstar,0,0.5,country);
	it_min=argmin(prix1+prix2);			# locate where the min of deaths+cost is
	copt=c[:,:,it_min];					# optimal control
	prixopt=Int(floor(prix2[it_min]));	# and its associative cost price (integer)

	if country=="France"
		S=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655];
	elseif country=="Burkina"
		S=18450494/100*[13.2,12.6,11.6,12.9,11.5,9.3,7.2,5.4,4.3,3.1,2.5,1.8,1.4,0.9,0.7,0.4,0.3,0.2];
	elseif country=="Vietnam"
		S=[4539031,7525542,6976064,6474991,7145151,8741274,8345522,7631483,6941367,6442420,5774836,5136526,4149565,2772903,1526235,1113819,873832,606471,268159,124992];
	end

	sumS=sum(S);
	B=Bstar*S/sumS;

	# The two strategies: (the cost of measure control being prixopt to distribute)  #

		# First strategy (youngs: less than 25 years old).

	duration=Int(floor(1+prixopt/(dt*sum(B[1:5])*(0.95^2)))); # the duration of the young control
	begincontrol=Int(1+floor(7/dt)) # start the control at the beginning of the second week (7 days)
	c1=zeros(nbr,Na);
	c1[begincontrol:(begincontrol+duration-1),1:5]=0.95*ones(duration,5);

	p1=plot(t,a*5,c1',st=:contourf,xlabel="Time (days)",ylabel="Age (years)",title="Control of <25s",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,clims=(0,1),xticks=collect(0:50:T),yticks=collect(0:20:100));

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,c1,Hsat,0.5,5.2,0.1,country);

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="Control of <25s",color=couleurs[1],width=2);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="Control of <25s",color=couleurs[1],width=2);

	TotalMort[2]=sum(Mdir[end,:])+sum(Mindir[end,:])
	ptot[2]=prop_totale[end];

	p5=plot!(p5,collect(0:5:((amax)*5)),[prop_age[end,:];prop_age[end,end]],label="Control of <25s",color=couleurs[1],width=2,linetype=:steppost);
	p5=plot!(p5,collect(0:5:(amax)*5),prop_totale[end]*ones(Na),label="",linestyle=:dot,color=couleurs[1],width=2);

		# Second strategy (age uniformly distributed). Let: 1+43491/(dt*sum(B)*(0.95^2))

	duration=Int(floor(1+prixopt/(dt*sum(B)*(0.95^2)))); # the duration of the uniform control
	c2=zeros(nbr,Na);
	c2[begincontrol:(begincontrol+duration-1),:]=0.95*ones(duration,Na);

	p2=plot(t,a*5,c2',st=:contourf,xlabel="Time (days)",ylabel="Age (years)",title="Uniform control",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,clims=(0,1),xticks=collect(0:50:T),yticks=collect(0:20:100));

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,c2,Hsat,0.5,5.2,0.1,country);

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="Uniform control",color=couleurs[2],width=2);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="Uniform control",color=couleurs[2],width=2);

	TotalMort[3]=sum(Mdir[end,:])+sum(Mindir[end,:])
	ptot[3]=prop_totale[end];

	p5=plot!(p5,collect(0:5:((amax)*5)),[prop_age[end,:];prop_age[end,end]],label="Uniform control",color=couleurs[2],width=2,linetype=:steppost);
	p5=plot!(p5,collect(0:5:(amax)*5),prop_totale[end]*ones(Na),label="",linestyle=:dot,color=couleurs[2],width=2);

		# Optimal plots

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,copt,Hsat,0.5,5.2,0.1,country);

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="Optimal control",color=couleurs[3],width=2);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="Optimal control",color=couleurs[3],width=2);
	p3=plot!(p3,t,log10.(tt+Hsat*ones(size(t))),linestyle=:dot,label="Healthcare capacity",color="red",legendfont = font("Times new roman", 24),legend=:right,width=2);

	TotalMort[4]=sum(Mdir[end,:])+sum(Mindir[end,:])
	ptot[4]=prop_totale[end];

	p5=plot!(p5,collect(0:5:((amax)*5)),[prop_age[end,:];prop_age[end,end]],label="Optimal control",color=couleurs[3],width=2,linetype=:steppost);
	p5=plot!(p5,collect(0:0.1:(amax)*5),(1-(1/3.3))*ones(Na),linestyle = :dash,label="Herd immunity",color="black",legend=:bottomleft,width=2);
	p5=plot!(p5,collect(0:5:(amax)*5),prop_totale[end]*ones(Na),label="",linestyle=:dot,color=couleurs[3],width=2);

	println("The function took $(time()-timecomput) seconds to execute.")

	return p1,p2,p3,p4,p5,TotalMort,ptot;

end

		##### Figures for practice measures control - Figure S4 - #####

function Practice(T)

	timecomput=time();	# initial time

	Bstar=1000;	TotalMort=zeros(3);	ptot=zeros(3);

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,0,5000,0.5,5.2,0.1,"France");

	dt=0.2;	nbr=Int(T/dt+1);
	Na=length(a);
	hosp_max=maximum(Hosp);
	couleurs=[:darkorange,:limegreen];

	tt=ones(length(t));
	aa=ones(Na);

	p3=plot(t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),xlabel="Time (days)",ylabel="\$\\log_{10}\$(Number of cases)",title="Number of hospitalisations",titlefontsize=38,xtickfont=font(30),ytickfont=font(30),xguidefontsize=32,yguidefontsize=30,xticks=collect(0:30:365),color="blue",label="Without control",legendfont = font("Times new roman", 26),width=2);

	p4=plot(a*5,log10.(aa+(Mindir)[nbr,:]+(Mdir)[nbr,:]),xlabel="Age (years)",ylabel="\$\\log_{10}\$(Deaths)",title="Cumulative number of deaths at t=$(T)",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,xticks=collect(0:20:100),color="blue",label="Without control",legendfont = font("Times new roman", 28),width=2);

	TotalMort[1]=sum(Mdir[end,:])+sum(Mindir[end,:]);
	ptot[1]=prop_totale[end];

		## Optimal strategy

	c,err,pp,prix,prix2=CovidOptim(210,1000,0,0.5,"France"); # the argmin is reached before 210th iterations so we stop there

	it_min=argmin(prix+prix2); # should be it_min=207;
	c_opt=c[:,:,it_min];

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,c_opt,5000,0.5,5.2,0.1,"France");

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="Optimal control",color=couleurs[1],width=2);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="Optimal control",color=couleurs[1],width=2);

	TotalMort[2]=sum(Mdir[end,:])+sum(Mindir[end,:])
	ptot[2]=prop_totale[end];

		## Discretised optimal strategy

	c_discret=c_opt;
	jours_total=21;						# window discretisation for time, in days
	age_total=2;						# window discretisation for age
	nbr_age=Int(Na/age_total);			# number of discretised age
	nbr_jours=Int(jours_total/dt);		# number of discretised days
	nbr_it=Int(floor(nbr/nbr_jours));	# number of discretised iterations
	discret_time=zeros(Na);
	discrete_age=zeros(Na);

	for ii=1:nbr_it

		discret_time=(sum(c_opt[(ii-1)*(nbr_jours)+1:ii*(nbr_jours),:],dims=1)[1,:])/(nbr_jours);

		for aak=1:nbr_age
			discrete_age[(aak-1)*age_total+1:aak*age_total]=(sum(discret_time[(aak-1)*age_total+1:aak*age_total])*ones(age_total))/(age_total);
		end

		c_discret[((ii-1)*(nbr_jours)+1):(ii*(nbr_jours)),:]=hcat(fill.(discrete_age,nbr_jours)...);
	end

		taille_disc=nbr-nbr_it*nbr_jours;

	for aak=1:nbr_age
		c_discret[nbr_it*nbr_jours+1:end,(aak-1)*age_total+1:aak*age_total]=(sum(c_opt[nbr_it*nbr_jours+1:end,(aak-1)*age_total+1:aak*age_total]))*ones(taille_disc,age_total)/(age_total*taille_disc);

	end

	pdisc=plot(t,a*5,c_discret',st=:contourf,xlabel="Time (days)",ylabel="Age (years)",title="Discretised control",titlefontsize=28,xtickfont=font(30),ytickfont=font(30),xguidefontsize=34,yguidefontsize=34,clims=(0,1),xticks=collect(0:50:T),yticks=collect(0:20:100));

	# Uniform cost for control measures

	SS=[3671719,4084036,4187992,4140996,3757482,3713426,4056469,4231788,4072226,4512223,4425730,4359376,4099662,3899944,3477098,2216562,1869006,1375537,678776,233655];
	sumSS=sum(SS);
	B=Bstar*SS/sumSS;
	prixx=dt*sum((c_discret.^2)*B);

	t,a,S,NewHosp,NewInfec,Hosp,NonHosp,Mdir,Mindir,R,prop_age,prop_totale=Covid(T,c_discret,5000,0.5,5.2,0.1,"France");

	p3=plot!(p3,t,log10.(tt+sum(Hosp,dims=2)[1:nbr,1]),label="Discretised control",color=couleurs[2],width=2,linestyle=:dot);
	p4=plot!(p4,a*5,log10.(aa+(Mdir)[nbr,:]+(Mindir)[nbr,:]),label="Discretised control",color=couleurs[2],width=2,linestyle=:dot);

	TotalMort[3]=sum(Mdir[end,:])+sum(Mindir[end,:])
	ptot[3]=prop_totale[end];
	p3=plot!(p3,t,log10.(tt+5000*ones(size(t))),linestyle=:dot,label="Healthcare capacity",color="red",legendfont = font("Times new roman", 24),legend=:topright,width=2);

	p5=plot(t/7,c_discret[:,1],xlabel="Time (weeks)",ylabel="Intensity of control",title="Discretised optimal control",titlefontsize=38,xtickfont=font(30),ytickfont=font(30),xguidefontsize=32,yguidefontsize=32,xticks=collect(0:3:55),yticks=collect(0:0.1:1),label="0-10 years",legendfont = font("Times new roman", 26),width=2);

	iii=1;

	while iii<6
		p5=plot!(p5,t/7,c_discret[:,2*iii+1],label="$(iii*10+1)-$((iii+1)*10) years",width=2);
		iii=iii+1;

	end

	p6=plot(t/7,c_discret[:,11],xlabel="Time (weeks)",ylabel="Intensity of control",title="Discretised optimal control",titlefontsize=38,xtickfont=font(30),ytickfont=font(30),xguidefontsize=32,yguidefontsize=32,xticks=collect(0:3:55),yticks=collect(0:0.1:1),legendfont = font("Times new roman", 26),width=2,label="$(iii*10+1)-$((iii+1)*10) years");
	iii=iii+1;

	while iii<10
		p6=plot!(p6,t/7,c_discret[:,2*iii+1],label="$(iii*10+1)-$((iii+1)*10) years",width=2);
		iii=iii+1;
	end

	println("The function took $(time()-timecomput) seconds to execute.")

	return pdisc,p3,p4,p5,p6,TotalMort,prixx,ptot;

end
