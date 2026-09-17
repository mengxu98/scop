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

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED };

struct cibersort_svm_parameter
{
	int svm_type;
	int kernel_type;
	int degree;
	double gamma;
	double coef0;

	double cache_size;
	double eps;
	double C;
	int nr_weight;
	int *weight_label;
	double* weight;
	double nu;
	double p;
	int shrinking;
	int probability;
};

struct cibersort_svm_model
{
	struct cibersort_svm_parameter param;
	int nr_class;
	int l;
	struct cibersort_svm_node **SV;
	double **sv_coef;
	double *rho;
	double *probA;
	double *probB;
	int *sv_indices;


	int *label;
	int *nSV;
	int free_sv;
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
void cibersort_svm_clear_diagnostics(void);
void cibersort_svm_flush_diagnostics(void);


#ifdef __cplusplus
}
#endif

void cibersort_svm_set_print_string_function(void (*print_func)(const char *));

#endif
