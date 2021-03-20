/*
forAll(designSpaceCells,i){
	const label j = designSpaceCells[i];
	scalar valEta = eta[j];
	//scalar val = alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty);
	alpha[j] = (alpha_s+(alpha_f-alpha_s)*valEta*(1.0+penalty)/(valEta + penalty));
	if(TempSetting) {
		DT[j] = (DTs+(DTf-DTs)*valEta*(1.0+penalty)/(valEta + penalty);
	}
}*/
