1.nucleotide  核苷酸变化图
  1.1 先将maf文件根据所需样本的分组，将各组的突变分别提取出来；
  1.2 然后利用Calculate.py去统计各核苷酸的变化，共6种；
      规则是：1.2.1 仅提取单核苷酸突变，不要indel；
	          1.2.3 分别计算6中核苷酸突变的比例：
			        I:   C->A and G->T
					II:  C->G and G->C
					III: C->T and G->A
					IV:  T->A and A->T
					V：  T->C and A->G
					VI： T->G and A->C
  1.3 然后，使用画百分比图，这个图可以直接用excel画出来；


2.trinucleotide  三联核苷酸变化图：
  2.1 先将maf文件根据所需样本的分组，将各组的突变分别提取出来；
  2.2 然后利用Calculate.py去统计各核苷酸的变化，共96种；
      具体有哪96种，可看“trinucleotide.xlsx”文件；
	  规则是：2.2.1 仅提取单核苷酸突变，不要indel；
	          2.2.2 分别将突变前后两个碱基补齐，然后再去统计各种碱基的比例：
			        可直接使用脚本“Calculate.py”
  2.3 使用“histogram_plot.R”作图，作图的input文件为上一步Calculate.py的输出文件；
      

3.Signature图
  3.1 用maf文件做signature的分析
  3.2 根据提供的分组信息，将每组样本的signature结果提取出来
  3.3 将signature结果进行作图，使用"hsit_plot_2.R"：
      格式如示例文件：TMB_upper_25.txt；
	  注意signature结果标题"Sample Name"在R中作图时，
	  空格应改为“SampleName”,否则会被识别成2个字符而报错;
	  
4. T.test分析
  4.1 用maf文件做signature的分析
  4.2 根据提供的分组信息，将每组样本的signature结果提取出来
  4.3 将signature结果进行T检验分析，使用脚本“T.test.R”
      具体每个signature改为什么title，脚本已转化，也可对应表：trinucleotide.xlsx
      
       