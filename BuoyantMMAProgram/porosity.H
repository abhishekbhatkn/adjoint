{
	forAll(eta,i){
		scalar valEta = eta[i];
		//scalar val = alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty);
		porosity[i] = (porosity_s+(porosity_f-porosity_s)*valEta*(1.0+penalty)/(valEta + penalty));
	}
	porosity.correctBoundaryConditions();
	
	//alphaTemp[j] = (porosity_s+(porosity_f-porosity_s)*valEta*(1.0+penalty)/(valEta + penalty));
	//DT = (DTs+(DTf-DTs)*eta*(1.0+penalty)/(eta + penalty));

}
