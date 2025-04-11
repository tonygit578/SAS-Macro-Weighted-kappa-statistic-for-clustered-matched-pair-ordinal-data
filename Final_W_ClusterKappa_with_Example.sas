/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
### SAS macro to calculate the weighted Kappa statistic and its variance under clustered data
 		Code to create the results in Examples section of the manuscript;
### "Weighted Kappa statistic for clustered matched-pair ordinal data";
### for Computational Statistics and Data Analysis;
### by Zhao Yang and Ming Zhou ;
### Version Date: 04Feb2014
** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;
/** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~** 
##### Revision History
On 24Sep2014
Based on the comments from Dr Sophie Vanbelle, update the formula by changing
2 # PV_ * W_M * PV_` * BigOmega * P_V * W_M * P_V`
to
2 # P_V * W_M * P_V` * BigOmega * PV_ * W_M * PV_`
in the component V22, however, there are no changes to the results

On 25Sep2014
More update on the the formula calculation
2 # P_V * W_M * P_V` * BigOmega * PV_ * W_M * PV_`
to
PV_ * W_M * P_V` * BigOmega * PV_ * W_M * P_V` + P_V * W_M * PV_` * BigOmega * P_V * W_M * PV_`
in the component V22
** ~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~~*~**/;

%macro W_ClusterKappa (dsin = , /* the input dataset */
					 clusV = , /* the cluster variable, need to be numeric*/
					 unitV = , /* the unit variable, need to be numeric*/
					 Resp1V =, /* Response variable for procedure 1, Row variable
					 			need to be numeric integer and > 0 */
					 Resp2V =,  /* Response variable for procedure 2, Column Variable
								need to be numeric integer and > 0 */ 
					 Weight = CA /* FC = Fleiss-Cohen weights; 
								CA = Cicchetti-Allison weights (default) */);
	proc sql noprint;
		create table clustDS as
			select distinct &clusV from &dsin
			order by &clusV; quit;
	data clustDS;
		set clustDS; clusterID + 1;
	proc sort data = &dsin;
		by &clusV &unitV;
	data resp_;
		merge &dsin clustDS;
		by &clusV;
	proc sort data = resp_;
		by clusterID &unitV;
	data F_Ana_DS (drop = &clusV &unitV rename = (clusterID = cluster));
		set resp_;
		by clusterID &unitV;
		if first.clusterID then subject = 0;
			subject + 1; 
	data anal_D;
		set F_Ana_DS (rename = (&Resp1V = resp1__ &Resp2V = resp2__ ));
		maxR = max(resp1__, resp2__);
	run;

	proc sql noprint;
		select max(maxR) into: N_Cat /*** to get the # of response category */ 
			from anal_D;
		select max(cluster) into: nclus /*** to get the # of clusters ***/
			from anal_D;
		create table count1 as  /** this is the calculate the marginal total for each cluster of treatment 1 */
			select distinct cluster, resp1__, count(*) as count1
				from anal_D
				group by cluster, resp1__;  
		create table count2 as  /** this is the calculate the marginal total for each cluster of treatment 2 */
			select distinct cluster, resp2__, count(*) as count2
				from anal_D
				group by cluster, resp2__; 
	quit; 
	
	** define the dataset containing weight;
	data weight;
		%do i = 1 %to &N_Cat;
			var&i = .;
		%end; output;
		data weight(drop = i j);
			set weight;
			array varlist{*} var1 -- var%sysfunc(left(&N_Cat));
			do i = 1 to &N_Cat;
				do j = 1 to &N_Cat;
					%if &Weight = FC %then %do;
						varlist{j} = 1 - (j-i)**2/(&N_Cat-1)**2;
					%end;
					%else %if &Weight = CA %then %do;
						varlist{j} = 1 - abs(j-i)/(&N_Cat-1) ;
					%end;
				end;
				output;
			end;
		run;		

	*** this is to construct the dataset with 0 in each combination;
	data constru1; /*** generate a dataset with full combination of tabulation */
		do cluster = 1 to &nclus;
			do resp1__ = 1 to &N_Cat;
				count1 = 0; output;
			end;
		end;
	proc sort data = constru1;  
		by cluster resp1__;
	data count1;	/** merge the datasets to get the full combination of tabulation */
		merge constru1 count1 ;
		by cluster resp1__;

	data constru2;
		do cluster = 1 to &nclus;
			do resp2__ = 1 to &N_Cat;
				count2 = 0; output;
			end;
		end;
	proc sort data = constru2;
		by cluster resp2__;
	data count2; /** merge the datasets to get the full combination of tabulation */
		merge constru2 count2 ;
		by cluster resp2__;
	run; quit;

	proc summary data = count1; /*** contains the # of total subjects */
		class cluster;
		var count1;
		output out = total (drop = _freq_ where = (_TYPE_ = 0) ) sum = total ;
	data _null_;
		set total (keep = total);
		call symput("total", total);
	run;

	proc summary data = count1; /*** contains the # of total subjects for each cluster;*/
		class cluster;
		var count1;
		output out = n_k (drop = _freq_ where = (_TYPE_ = 1) ) sum = n_k ;
	proc sort data = n_k (drop = _type_) out = n_k;
		by cluster;
	run;  

	/** to calculate the count for each cell in the kxk table;*/
	ods listing close;
		proc freq data = anal_D;
			by cluster; 
			tables resp1__ * resp2__/out = pctc (drop = PERCENT);
		run;
	ods listing;
	proc sort data = pctc;
		by cluster resp1__ resp2__;
	run;

	data fake_D;
		do cluster = 1 to &nclus;
			do resp1__ = 1 to &N_Cat;
				do resp2__ = 1 to &N_Cat;
					count = 0;
					output;
				end;
			end;
		end;
	proc sort data = fake_D;
		by cluster resp1__ resp2__;
	run;
			
	data pctc;
		merge fake_D pctc;
		by cluster resp1__ resp2__;
	run;
	data pctc;
		set pctc;
		by cluster resp1__ resp2__;
		if first.cluster then clseq = 0;
			clseq + 1;
	proc sql noprint;
		select max(clseq) into: maxSeq
		from pctc; quit;
	proc sort data = pctc;
		by cluster clseq;
	proc transpose data =  pctc out = concord (drop = _LABEL_ _NAME_) prefix = nn;
		by cluster;
		id clseq;
		var count;
	proc sort data = concord;
		by cluster;
	run;  

	%macro outputM(countds = );	
		%do i = 1 %to &N_Cat;
			%let seqD = %sysfunc( substr( &countds, 6, 1));
			data margin&seqD.&i (rename = (count&seqD = %if &seqD = 1 %then %do; n&i._k %end; %else %do; n_&i.k %end; ) drop = resp&seqD.__); 
				set &countds;
				if resp&seqD.__ = &i then output;
			proc sort data = margin&seqD.&i;
				by cluster;
			run;
		%end;
	%mend;
	%outputM(countds = Count1); /*** output the marginal total for treatment 1*/;
	%outputM(countds = Count2); /*** output the marginal total for treatment 2*/;

	data finalD;
		merge N_k 
			%do i = 1 %to &N_Cat;
				margin1&i margin2&i 
			%end;	concord;
		by cluster;
	data finalD;
		set finalD;
		BigN = &total;
	run;

	data finalD_A (drop = i);
		set finalD;
		array varlistA {*} %do i = 1 %to &N_Cat; n&i._k n_&i.k %end; nn1 - nn%sysfunc(left(&maxSeq));
		array varlistB {*} %do i = 1 %to &N_Cat; p&i._k p_&i.k %end; pp1 - pp%sysfunc(left(&maxSeq));
		do i = 1 to dim(varlistA);
			varlistB {i} = varlistA {i}/n_k;
		end;
		omega_k = n_k/BigN;
	run;
			
	*** calculation ussing the Kappa under the independence (non-clusterd situation);
	data fakestr;
		do resp1__ = 1 to &N_Cat;
			do resp2__ = 1 to &N_Cat;
				count = 0;	output;
			end;
		end;
	proc sort data = fakestr;
		by resp1__ resp2__;
	run;

	proc summary data = anal_D;
		class resp1__ resp2__;
		output out = interm_D (where = (_TYPE_ = 3) rename = (_FREQ_ = count) );
	proc sort data = interm_D (drop = _TYPE_) out = interm_D;
		by resp1__ resp2__;
	data A_Response;
		merge fakestr interm_D;
		by resp1__ resp2__;
	run;

	title2 "*** Part 1: weighted Kappa statistics assuming independence *** ";	
	proc freq data = A_Response;
		tables resp1__ * resp2__/%if &Weight = FC %then %do;
									agree(WT = &Weight);
								 %end; %if &Weight = CA %then %do;
									agree ;
								 %end;
		weight count/zeros;
	run;

	proc iml;
		use weight; read all ;
		W_M = var1 %do i = 2 %to &N_Cat; || var&i %end; ;
		close weight;
		* Read data into IML ;
		use finalD_A;  read all ;

		kappOD = J( 1, 6, .);
		*** this is used to calculate the point estimation;
		%do i = 1 %to &N_Cat;
			  p&i._ = omega_k` * p&i._k;
			  p_&i = omega_k` * p_&i.k;
		%end;
		%do i = 1 %to &maxSeq;
			  pp&i._ = omega_k` * pp&i;
		%end;
		%do i = 1 %to &N_Cat;
			  %let start = %eval( (&i - 1) * &N_Cat + 1);
			  %let start1 = %eval(&start + 1);
			  %let start2 = %eval(&i * &N_Cat);
			  p_Part&i = pp&start._ %do j = &start1 %to &start2; || pp&j._ %end; ; 
		%end;

		p_M = p_Part1 %do i = 2 %to &N_Cat; // p_Part&i %end; ;
		pi_ = P1_ %do i = 2 %to &N_Cat; || P&i._ %end; ;
		p_i = P_1 %do i = 2 %to &N_Cat; || P_&i %end; ;

		Prod_WP = W_M # p_M; 
		P_O = sum(Prod_WP);
		P_E = pi_ * W_M * t(p_i); 

		kappa = (P_O - P_E) / (1 - P_E); ** this is the calculated Kappa value;

		*** end of point estimation calculation;
		*** this part is for the calculation of the variance based on Delta method;
		PO_k = pp1 %do i = 2 %to &maxSeq; ||pp&i %end; ;

		W_M_H = W_M[1,] %do i = 2 %to &N_Cat; ||W_M[&i,] %end; ;
		PO_V = PO_k * W_M_H`;

		P_V  = p1_k %do i = 2 %to &N_Cat; || p&i._k %end; ;
		PV_  = p_1k %do i = 2 %to &N_Cat; || p_&i.k %end; ;

		omega_k2 = omega_k # omega_k;
		BigOmega = ( I(&nclus) - omega_k * J(1, &nclus, 1) ) * diag(omega_k2) * ( I(&nclus) - J(&nclus, 1, 1) * t(omega_k) );
			  
		V11 = PO_V` * BigOmega * PO_V;
		V12 = PO_V` * BigOmega * (P_V * W_M * PV_` + PV_ * W_M * P_V`) * omega_k;
		V21 = V12;
		V22 = omega_k` * (PV_ * W_M * P_V` * BigOmega * P_V * W_M * PV_` + P_V * W_M * PV_` * BigOmega * PV_ * W_M * P_V` 
							+ PV_ * W_M * P_V` * BigOmega * PV_ * W_M * P_V` + P_V * W_M * PV_` * BigOmega * P_V * W_M * PV_`) * omega_k;

		VMatrix = J( 2, 2, .);
		LVector = J( 1, 2, .);

		VMatrix[1,1] = V11;
		VMatrix[1,2] = V12;
		VMatrix[2,1] = V21;
		VMatrix[2,2] = V22;
		VMatrix = (&nclus # &nclus) # VMatrix/(&nclus - 1);
			  
		LVector[1] = 1 / (1 - P_E);
		LVector[2] = (P_O - 1) / ( (1 - P_E) # (1 - P_E) );

		ASE_Del = sqrt( LVector * VMatrix * LVector`/&nclus ); *** variance using Delta method;
			  
		lower = kappa - 1.96 # ASE_Del;
		upper = kappa + 1.96 # ASE_Del;
			 * print V11, V12, V21, V22 VMatrix Var_Del P_O kappa;

			  *** end of the calculation of the variance based on Delta method;

			  ******* this part is for pooling the final results information;
		
		title2 "*** Part 2: Weighted Kappa statistics for clustered matched-pair ordinal data *** ";	
		title3 "************************  Using Weighting Scheming: &weight ********************* ";	
		print P_O P_E kappa ASE_Del[label="Standard Error"],, 
			  lower[label="95% CI Lower limit"] upper[label="95% CI Upper limit"];
		quit;
%mend;



/***** Here is the example ********/;
data cerv;
  input id y11 y12 y21 y22;
  cards;
1 4 4 3 4
2 1 1 1 1
3 2 1 1 1
4 3 3 2 3
5 1 2 1 3
6 1 1 1 1
7 2 2 2 2
8 3 2 2 2
9 4 3 3 2
10 2 2 1 1
11 2 1 1 1
12 4 4 3 2
13 4 2 1 1
14 1 2 1 1
15 4 4 4 4
16 4 4 4 4
17 4 2 2 2
18 2 2 1 1
19 2 2 1 1
20 4 2 2 2
21 4 4 4 3
22 4 3 3 3
23 1 1 1 1
24 4 1 1 1
25 2 1 2 1
26 4 2 2 2
27 4 4 4 4
28 2 2 1 1
29 2 2 3 1
30 4 3 4 4
31 4 4 4 4
32 2 3 2 2
33 3 1 2 2
34 4 3 2 2
35 1 1 1 1
36 1 1 1 1
37 4 3 2 2
38 1 1 1 1
39 1 1 2 2
40 4 4 4 4
41 3 2 2 2
42 4 3 2 2
43 1 1 1 1
44 1 1 1 1
45 4 4 4 4
46 2 2 1 1
47 3 1 2 3
48 2 2 1 1
49 3 2 2 2
50 3 2 3 2
51 4 3 4 3
52 1 1 1 1
53 2 2 2 1
54 2 2 2 1
55 2 1 1 1
56 4 3 2 2
57 4 4 4 4
58 4 3 3 2
59 3 3 2 2
60 4 3 1 1
61 2 2 2 1
62 3 2 2 1
63 2 1 1 1
64 2 1 1 1
65 2 2 1 1
66 4 3 2 1
67 1 1 2 2
68 1 1 2 1
69 2 2 1 1
70 2 1 1 1
71 2 3 2 2
72 2 1 2 2
73 1 1 1 1
74 2 2 2 2
75 4 4 4 2
76 3 2 2 2
77 4 3 2 2
78 2 3 2 2
79 3 1 2 3
80 2 1 1 1
81 3 3 2 2
82 2 1 1 2
83 2 2 1 1
84 3 2 2 2
85 2 2 2 2
;
run;
data cerv;
	set cerv;
	trt1 = 1; trt2 = 2;
run;
data response;
	set cerv (rename = (id = cluster));
	array ary1{*} y11 y12;
	array ary2{*} y21 y22;
	do i = 1 to dim(ary1);
		Resp1 =  ary1{i} ;
		Resp2 =  ary2{i} ;
		Subject = i;
		output;
	end;
run;

%W_ClusterKappa(dsin = response,  clusV = cluster,  unitV = subject,  Resp1V = Resp1, Resp2V = Resp2, weight = FC);
%W_ClusterKappa(dsin = response,  clusV = cluster,  unitV = subject,  Resp1V = Resp1, Resp2V = Resp2, weight = CA);

