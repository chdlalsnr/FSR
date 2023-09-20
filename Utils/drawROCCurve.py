from sklearn.metrics import roc_curve
y_pred = model.predict(x_test).ravel()
fpr,tpr,threshold = roc_curve(y_test, y_pred)

from sklearn.metrics import auc
auc = auc(fpr, tpr)

# Plot ROC curve
plt.plot([0, 1],[0, 1], 'k--',alpha=0.6)
plt.plot(fpr, tpr, c="red", label='AUC = {:.3f})'.format(auc))

plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('HH Classifier Model ROC Curve')
plt.legend(loc='best')
plt.grid(True, axis='y', alpha=0.5, linestyle='--')

plt.show()

# Draw accuracy and validation loss
y_vloss = history.history['val_loss']
y_acc  = history.history['acc']

x_len = np.arange(len(y_acc))
plt.plot(x_len, y_vloss, "o", c="red",  markersize=3, label="Validation Loss")
plt.plot(x_len, y_vloss, c="red", alpha=0.5, linewidth=1)
plt.plot(x_len, y_acc  , "o", c="blue", markersize=3, label="Accuracy")
plt.plot(x_len, y_acc, c="blue", alpha=0.5, linewidth=1)
plt.legend()

plt.title("Test Model Accuracy and Validation Loss")
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.grid(True, axis='y', alpha=0.5, linestyle='--')

plt.show()
