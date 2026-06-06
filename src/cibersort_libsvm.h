#ifndef _SCOP_CIBERSORT_LIBSVM_H
#define _SCOP_CIBERSORT_LIBSVM_H

#define LIBSVM_VERSION 323

#ifdef __cplusplus
extern "C" {
#endif

extern int cibersort_libsvm_version;

struct cibersort_svm_node
{
	int index;
	double value;
};

struct cibersort_svm_problem
{
	int l;
	double *y;
	struct cibersort_svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct cibersort_svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;	/* for poly */
	double gamma;	/* for poly/rbf/sigmoid */
	double coef0;	/* for poly/sigmoid */

	/* these are for training only */
	double cache_size; /* in MB */
	double eps;	/* stopping criteria */
	double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
	int nr_weight;		/* for C_SVC */
	int *weight_label;	/* for C_SVC */
	double* weight;		/* for C_SVC */
	double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
	double p;	/* for EPSILON_SVR */
	int shrinking;	/* use the shrinking heuristics */
	int probability; /* do probability estimates */
};

/*
//
// cibersort_svm_model
//
*/
struct cibersort_svm_model
{
	struct cibersort_svm_parameter param;	/* parameter */
	int nr_class;		/* number of classes, = 2 in regression/one class svm */
	int l;			/* total #SV */
	struct cibersort_svm_node **SV;		/* SVs (SV[l]) */
	double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
	double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
	double *probA;		/* pariwise probability information */
	double *probB;
	int *sv_indices;        /* sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */

	/* for classification only */

	int *label;		/* label of each class (label[k]) */
	int *nSV;		/* number of SVs for each class (nSV[k]) */
				/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
	/* XXX */
	int free_sv;		/* 1 if cibersort_svm_model is created by cibersort_svm_load_model*/
				/* 0 if cibersort_svm_model is created by cibersort_svm_train */
};

struct cibersort_svm_model *cibersort_svm_train(const struct cibersort_svm_problem *prob, const struct cibersort_svm_parameter *param);
void cibersort_svm_cross_validation(const struct cibersort_svm_problem *prob, const struct cibersort_svm_parameter *param, int nr_fold, double *target);

int cibersort_svm_save_model(const char *model_file_name, const struct cibersort_svm_model *model);
struct cibersort_svm_model *cibersort_svm_load_model(const char *model_file_name);

int cibersort_svm_get_svm_type(const struct cibersort_svm_model *model);
int cibersort_svm_get_nr_class(const struct cibersort_svm_model *model);
void cibersort_svm_get_labels(const struct cibersort_svm_model *model, int *label);
void cibersort_svm_get_sv_indices(const struct cibersort_svm_model *model, int *sv_indices);
int cibersort_svm_get_nr_sv(const struct cibersort_svm_model *model);
double cibersort_svm_get_svr_probability(const struct cibersort_svm_model *model);

double cibersort_svm_predict_values(const struct cibersort_svm_model *model, const struct cibersort_svm_node *x, double* dec_values);
double cibersort_svm_predict(const struct cibersort_svm_model *model, const struct cibersort_svm_node *x);
double cibersort_svm_predict_probability(const struct cibersort_svm_model *model, const struct cibersort_svm_node *x, double* prob_estimates);

void cibersort_svm_free_model_content(struct cibersort_svm_model *model_ptr);
void cibersort_svm_free_and_destroy_model(struct cibersort_svm_model **model_ptr_ptr);
void cibersort_svm_destroy_param(struct cibersort_svm_parameter *param);

const char *cibersort_svm_check_parameter(const struct cibersort_svm_problem *prob, const struct cibersort_svm_parameter *param);
int cibersort_svm_check_probability_model(const struct cibersort_svm_model *model);

//void cibersort_svm_set_print_string_function(void (*print_func)(const char *));

#ifdef __cplusplus
}
#endif

void cibersort_svm_set_print_string_function(void (*print_func)(const char *));

#endif /* _SCOP_CIBERSORT_LIBSVM_H */
