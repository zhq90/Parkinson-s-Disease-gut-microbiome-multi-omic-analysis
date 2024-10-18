Four ML algorithms were used

1. Logistic Regression

```python
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, RocCurveDisplay, f1_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
from sklearn import metrics

df = pd.read_csv("lassonet_top25_matrix.csv", header = 0, sep = ",")
md = pd.read_csv("metadataPDHP.csv", header = 0, sep = ",")
df = df.set_index("ID")
df = df.T
md = md.set_index("sampleID")
md = md.loc[df.index,"Group"]

lr = LogisticRegression(random_state=1,solver="saga", penalty = "l2",max_iter=100)
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=100, random_state = 99)

xiao_df = pd.read_csv("xiao_validation_matrix.csv", header = 0, sep = ",")
xiao_md = pd.read_csv("xiao_metadata.csv", header = 0, sep = ",")

xiao_df = xiao_df.set_index("Name")
xiao_df = xiao_df.loc[df.columns,]
xiao_md = xiao_md.set_index("Run")
xiao_md = xiao_md.loc[xiao_df.columns,:]
xiao_md = xiao_md["Group"]
xiao_df = xiao_df.T

tprs = []
fprs = []
aucs = []
mean_fpr = np.linspace(0,1,10000)
X = np.array(df)
y = np.array(md)
for i, (train, test) in enumerate(cv.split(X, y)):
    lr.fit(X[train],y[train])
    viz = RocCurveDisplay.from_estimator(
        lr,
        X[test],
        y[test],
        name = "Random Forest",
        alpha = 0.5,
        lw = 1,
        ax = None,
        pos_label = "PD",
        )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    fprs.append(np.mean(viz.fpr))


lr = LogisticRegression(random_state=1,solver="saga", penalty = "l2",max_iter=100)
lr.fit(df, md)
lr.score(xiao_df, xiao_md)
lr.score(df, md)

xiao_probs = lr.predict_proba(xiao_df)
probs = xiao_probs[:,1]
auc_value = roc_auc_score(xiao_md, probs)
auc_value = round(auc_value, 2)
auc_value
fpr, tpr, thresholds = roc_curve(xiao_md, probs, pos_label="PD")

validation_specificity = 1 - np.mean(fpr)
validation_sensitivity = np.mean(tpr)

fig, ax = plt.subplots()
ax.plot([0,1], [0,1], linestyle = "--", lw = 1, color = "#3A3A3A", alpha = 0.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="#37b7fa",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

ax.plot(
        fpr,
        tpr,
        color = "#fa7a37",
        label = r"Validation Set ROC (AUC = %0.2f)" % (auc_value),
        lw = 2,
        alpha = 0.8)

ax.set(xlim=[-0.05, 1.05],
       ylim=[-0.05, 1.05],
       title = "ROC curve (Logistic Regression)")
ax.legend(loc = "lower right", prop = { "size": 10 })
ax.set_xlabel("1 - Specificity")
ax.set_ylabel("Sensitivity")
plt.savefig("LogisticRegression_ROC_curve.pdf")
```

2. Random Forest

```python
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, RocCurveDisplay, f1_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
from sklearn import metrics

df = pd.read_csv("lassonet_top25_matrix.csv", header = 0, sep = ",")
md = pd.read_csv("metadataPDHP.csv", header = 0, sep = ",")
df = df.set_index("ID")
df = df.T
md = md.set_index("sampleID")
md = md.loc[df.index,"Group"]

rf = RandomForestClassifier(random_state=1, n_estimators = 150)
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=100, random_state = 99)

xiao_df = pd.read_csv("xiao_validation_matrix.csv", header = 0, sep = ",")
xiao_md = pd.read_csv("xiao_metadata.csv", header = 0, sep = ",")

xiao_df = xiao_df.set_index("Name")
xiao_df = xiao_df.loc[df.columns,]
xiao_md = xiao_md.set_index("Run")
xiao_md = xiao_md.loc[xiao_df.columns,:]
xiao_md = xiao_md["Group"]
xiao_df = xiao_df.T

tprs = []
fprs = []
aucs = []
mean_fpr = np.linspace(0,1,10000)
X = np.array(df)
y = np.array(md)
for i, (train, test) in enumerate(cv.split(X, y)):
    rf.fit(X[train],y[train])
    viz = RocCurveDisplay.from_estimator(
        rf,
        X[test],
        y[test],
        name = "Random Forest",
        alpha = 0.5,
        lw = 1,
        ax = None,
        pos_label = "PD",
        )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    fprs.append(np.mean(viz.fpr))


rf = RandomForestClassifier(n_estimators = 150, random_state=1)
rf.fit(df, md)
rf.score(xiao_df, xiao_md)
rf.score(df, md)

xiao_probs = rf.predict_proba(xiao_df)
probs = xiao_probs[:,1]
auc_value = roc_auc_score(xiao_md, probs)
auc_value = round(auc_value, 2)
auc_value
fpr, tpr, thresholds = roc_curve(xiao_md, probs, pos_label="PD")

validation_specificity = 1 - np.mean(fpr)
validation_sensitivity = np.mean(tpr)

fig, ax = plt.subplots()
ax.plot([0,1], [0,1], linestyle = "--", lw = 1, color = "#3A3A3A", alpha = 0.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="#37b7fa",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

ax.plot(
        fpr,
        tpr,
        color = "#fa7a37",
        label = r"Validation Set ROC (AUC = %0.2f)" % (auc_value),
        lw = 2,
        alpha = 0.8)

ax.set(xlim=[-0.05, 1.05],
       ylim=[-0.05, 1.05],
       title = "ROC curve (Random Forest)")
ax.legend(loc = "lower right", prop = { "size": 10 })
ax.set_xlabel("1 - Specificity")
ax.set_ylabel("Sensitivity")
plt.savefig("RandomForest_ROC_curve.pdf")
```

3. Neural Network

```python
import pandas as pd
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, RocCurveDisplay, f1_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
from sklearn import metrics

df = pd.read_csv("lassonet_top25_matrix.csv", header = 0, sep = ",")
md = pd.read_csv("metadataPDHP.csv", header = 0, sep = ",")
df = df.set_index("ID")
df = df.T
md = md.set_index("sampleID")
md = md.loc[df.index,"Group"]

nn = MLPClassifier(random_state=9, activation = "relu", max_iter = 2000, hidden_layer_sizes=(10000,))
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=100, random_state = 99)

xiao_df = pd.read_csv("xiao_validation_matrix.csv", header = 0, sep = ",")
xiao_md = pd.read_csv("xiao_metadata.csv", header = 0, sep = ",")

xiao_df = xiao_df.set_index("Name")
xiao_df = xiao_df.loc[df.columns,]
xiao_md = xiao_md.set_index("Run")
xiao_md = xiao_md.loc[xiao_df.columns,:]
xiao_md = xiao_md["Group"]
xiao_df = xiao_df.T

tprs = []
fprs = []
aucs = []
mean_fpr = np.linspace(0,1,10000)
X = np.array(df)
y = np.array(md)
for i, (train, test) in enumerate(cv.split(X, y)):
    nn.fit(X[train],y[train])
    viz = RocCurveDisplay.from_estimator(
        nn,
        X[test],
        y[test],
        name = "Random Forest",
        alpha = 0.5,
        lw = 1,
        ax = None,
        pos_label = "PD",
        )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    fprs.append(np.mean(viz.fpr))


nn = MLPClassifier(random_state=9, activation = "relu", max_iter = 2000, hidden_layer_sizes=(10000,))
nn.fit(df, md)
nn.score(xiao_df, xiao_md)
nn.score(df, md)

xiao_probs = nn.predict_proba(xiao_df)
probs = xiao_probs[:,1]
auc_value = roc_auc_score(xiao_md, probs)
auc_value = round(auc_value, 2)
auc_value
fpr, tpr, thresholds = roc_curve(xiao_md, probs, pos_label="PD")

validation_specificity = 1 - np.mean(fpr)
validation_sensitivity = np.mean(tpr)

fig, ax = plt.subplots()
ax.plot([0,1], [0,1], linestyle = "--", lw = 1, color = "#3A3A3A", alpha = 0.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="#37b7fa",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

ax.plot(
        fpr,
        tpr,
        color = "#fa7a37",
        label = r"Validation Set ROC (AUC = %0.2f)" % (auc_value),
        lw = 2,
        alpha = 0.8)

ax.set(xlim=[-0.05, 1.05],
       ylim=[-0.05, 1.05],
       title = "ROC curve (Neural Network)")
ax.legend(loc = "lower right", prop = { "size": 10 })
ax.set_xlabel("1 - Specificity")
ax.set_ylabel("Sensitivity")
plt.savefig("neuralnetwork_ROC_curve.pdf")
```

4. SVM

```python
import pandas as pd
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, auc, RocCurveDisplay, f1_score, classification_report, confusion_matrix
import matplotlib.pyplot as plt
from sklearn import metrics

df = pd.read_csv("lassonet_top25_matrix.csv", header = 0, sep = ",")
md = pd.read_csv("metadataPDHP.csv", header = 0, sep = ",")
df = df.set_index("ID")
df = df.T
md = md.set_index("sampleID")
md = md.loc[df.index,"Group"]

svc = SVC(kernel = "rbf", random_state = 9, probability = True)
cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=100, random_state = 99)

xiao_df = pd.read_csv("xiao_validation_matrix.csv", header = 0, sep = ",")
xiao_md = pd.read_csv("xiao_metadata.csv", header = 0, sep = ",")

xiao_df = xiao_df.set_index("Name")
xiao_df = xiao_df.loc[df.columns,]
xiao_md = xiao_md.set_index("Run")
xiao_md = xiao_md.loc[xiao_df.columns,:]
xiao_md = xiao_md["Group"]
xiao_df = xiao_df.T

tprs = []
fprs = []
aucs = []
mean_fpr = np.linspace(0,1,10000)
X = np.array(df)
y = np.array(md)
for i, (train, test) in enumerate(cv.split(X, y)):
    svc.fit(X[train],y[train])
    viz = RocCurveDisplay.from_estimator(
        svc,
        X[test],
        y[test],
        name = "Random Forest",
        alpha = 0.5,
        lw = 1,
        ax = None,
        pos_label = "PD",
        )
    interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    interp_tpr[0] = 0.0
    tprs.append(interp_tpr)
    aucs.append(viz.roc_auc)
    fprs.append(np.mean(viz.fpr))


svc = SVC(kernel = "rbf", random_state = 9, probability = True)
svc.fit(df, md)
svc.score(xiao_df, xiao_md)
svc.score(df, md)

xiao_probs = svc.predict_proba(xiao_df)
probs = xiao_probs[:,1]
auc_value = roc_auc_score(xiao_md, probs)
auc_value = round(auc_value, 2)
auc_value
fpr, tpr, thresholds = roc_curve(xiao_md, probs, pos_label="PD")

validation_specificity = 1 - np.mean(fpr)
validation_sensitivity = np.mean(tpr)

fig, ax = plt.subplots()
ax.plot([0,1], [0,1], linestyle = "--", lw = 1, color = "#3A3A3A", alpha = 0.8)

mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="#37b7fa",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

ax.plot(
        fpr,
        tpr,
        color = "#fa7a37",
        label = r"Validation Set ROC (AUC = %0.2f)" % (auc_value),
        lw = 2,
        alpha = 0.8)

ax.set(xlim=[-0.05, 1.05],
       ylim=[-0.05, 1.05],
       title = "ROC curve (SVM)")
ax.legend(loc = "lower right", prop = { "size": 10 })
ax.set_xlabel("1 - Specificity")
ax.set_ylabel("Sensitivity")
plt.savefig("SVM_ROC_curve.pdf")
```
