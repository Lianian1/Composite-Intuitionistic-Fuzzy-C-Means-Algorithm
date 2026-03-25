function EVAL = EvaluateFast(ACTUAL,PREDICTED)
PREDICTED = (calculate_true_labels(PREDICTED',ACTUAL))';

% confusion.getMatrix 第三个参数 Display=0
[~,Result,~]= confusion.getMatrix(ACTUAL',PREDICTED',0);

nmi = fNMI(PREDICTED',ACTUAL');
EVAL = [Result.Accuracy, Result.F1_score, nmi];
end
