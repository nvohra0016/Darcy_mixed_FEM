// Script to solve for Darcy flow using mixed finite elements.
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_raviart_thomas.h> // Raviart-Thomas fe is declared in this
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bdm.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/convergence_table.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/packaged_operation.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/base/timer.h>


#include <deal.II/fe/mapping_q_eulerian.h>

// Sparse direct solve
#include <deal.II/lac/sparse_direct.h>

// Sparse ILU
#include <deal.II/lac/sparse_ilu.h>


#include <fstream>
#include <iostream>

#include <math.h>

using namespace dealii;

class Darcyeq{
    public:
       Darcyeq();
       void run();

    private:
        void make_grid();
        void setup_system();
        void assemble_system();
        void solve();
        void output_results();

        Triangulation<2> triangulation;
        DoFHandler<2> dof_handler;
        FESystem<2> fe;
        unsigned int np, nqf; // DoFs of pressure p and flux qf

        BlockSparsityPattern sparsity_pattern, preconditioner_sparsity_pattern;
        BlockSparseMatrix<double> system_matrix, preconditioner_matrix;
        BlockVector<double> solution, system_rhs;

        AffineConstraints<double> constraints;

        // Grid discretization parameters
        double hx, hy;
};

Darcyeq::Darcyeq() : dof_handler(triangulation), fe(FE_RaviartThomas<2>(0), 1, FE_DGQ<2>(0), 1),
    hx(1.0/200), hy(1.0/200)
{}

// RHS external source function \nabla \cdot qf = f
class rhs_f: public Function<2> {
    public:
    
    rhs_f(): Function<2>(1)
    {}

    virtual double value (const Point<2> &p, const unsigned int component = 0) const override;
};

double rhs_f::value (const Point<2> &p, const unsigned int component) const {
    double t = this->get_time();
    double x = p[0], y = p[1];
    
    double val = 0.0;

    return val;
}


// Permeability 
class K_inverse: public TensorFunction<2,2> {
    public:
        K_inverse(): TensorFunction<2,2>()
        {}

       virtual void value_list(const std::vector<Point<2>> &points, std::vector<Tensor<2,2>> &values) const;
};

void K_inverse::value_list(const std::vector<Point<2>> &points, std::vector<Tensor<2,2>> &values) const {

    double val;
    double val1 = 1.e-6, val2 = 1.e-2;

    for(unsigned int p = 0; p < points.size(); ++p) {
        double x = points[p][0], y = points[p][1];
        
        values[p].clear();

        if (std::abs(x - 0.5) < 0.2 + 1.e-6 && std::abs(y - 0.5) < 0.2 + 1.e-6) {
            val = 1.0/val1;
        }
        else {
            val = 1.0/val2;
        }

        for(unsigned int d = 0; d < 2; ++d) {
            values[p][d][d] = val;
        }
    }
}

// Boundary values function
class pressure_boundary: public Function<2> {
    public: 

    pressure_boundary(): Function<2>(1)
    {}

    virtual double value (const Point<2> &p, const unsigned int component = 0) const override;
};

double pressure_boundary::value (const Point<2> &p, const unsigned int component) const {
    double t = this->get_time();
    double x = p[0], y = p[1];

    double val = 0.0;

    if (std::abs(x - 0.0) < 1.e-6) {
        val = 1.e6;
    } 

    return val;
}

// Flux boundary values
class flux_boundary : public Function<2> // Flux boundary values
{
public:
  flux_boundary() : Function<2>(2)
  {}
  
  virtual void vector_value (const Point<2> &p, Vector<double> &value) const;
};

void flux_boundary::vector_value (const Point<2> &p, Vector<double> &value) const
{
  double t = this->get_time();
  double x = p[0], y = p[1];
  double pi = numbers::PI;
  
  value(0) = 0.0;
  value(1) = 0.0;
}

class boundary_values : public Function<2> // Boundary values for displacements.
{
public:
  boundary_values () : Function<2>(5)
  {}
  
  virtual void vector_value (const Point<2> &p, Vector<double> &value) const;
};

void boundary_values::vector_value (const Point<2> &p, Vector<double> &value) const
{
  double t = this->get_time();
  double x = p[0], y = p[1];
  double pi = numbers::PI;

  double val_qx, val_qy, val_p;

  val_qx = 0.0;
  val_qy = 0.0;
  val_p = 0.0;

  value(0) = val_qx;
  value(1) = val_qy;
  value(2) = val_p;
}


//
// Main functions *--------------------------------*
//

// Generate grid
void Darcyeq::make_grid() {

    // Define grid cells and end points
    Point<2> p1, p2;
    p1[0] = 0.0; p1[1] = 0.0;
    p2[0] = 1.0; p2[1] = 1.0;

    std::vector<unsigned int> repititions(2);
    repititions[0] = static_cast<unsigned int>((p2[0] - p1[0])/hx);
    repititions[1] = static_cast<unsigned int>((p2[1] - p1[1])/hy);

    // Generate grid
    GridGenerator::subdivided_hyper_rectangle(triangulation, repititions, p1, p2);

    // Assign boundary labels
    for (const auto &cell: triangulation.cell_iterators()) {
        for (unsigned int face_n = 0; face_n < GeometryInfo<2>::faces_per_cell; ++face_n) {
            const auto center = cell->face(face_n)->center();

            if (std::fabs(center(1) - 0.0) < 1.e-6) {
                cell->face(face_n)->set_boundary_id(0);
            }
            else if (std::fabs(center(0) - 1.0) < 1.e-6) { 
                cell->face(face_n)->set_boundary_id(1);
            }
            else if (std::fabs(center(1) - 1.0) < 1.e-6) { 
                cell->face(face_n)->set_boundary_id(2);
            }
            else if (std::fabs(center(0) - 0.0) < 1.e-6) {
                cell->face(face_n)->set_boundary_id(3);
            }
        }
    }

    // Print info
    std::cout<<"Number of active cells =  "<<triangulation.n_active_cells()<<std::endl;
    std::cout<<"make_grid() completed"<<std::endl;
}

// Setup system
void Darcyeq::setup_system() {
    
    dof_handler.distribute_dofs(fe);
    DoFRenumbering::component_wise(dof_handler);

    const std::vector<types::global_dof_index> dofs_per_component = DoFTools::count_dofs_per_fe_component(dof_handler);
    nqf = dofs_per_component[0];
    np = dofs_per_component[2];

    std::cout<<"DoFs: Pressure np = "<<np<<"; Flux nqf = "<<nqf<<std::endl;

    const FEValuesExtractors::Scalar p(2); // dimension = 2
    const FEValuesExtractors::Vector qf(0);
    const FEValuesExtractors::Scalar qfx(0);
    const FEValuesExtractors::Scalar qfy(1);

    // Boundary constraints
    constraints.clear();

    flux_boundary flux_boundary_values_object;

    VectorTools::project_boundary_values_div_conforming(dof_handler, 0, flux_boundary_values_object, 0, constraints);
    //VectorTools::project_boundary_values_div_conforming(dof_handler, 0, flux_boundary_values_object, 1, constraints);
    VectorTools::project_boundary_values_div_conforming(dof_handler, 0, flux_boundary_values_object, 2, constraints);

    constraints.close();

    std::cout<<"Boundary constraints applied"<<std::endl;

    const std::vector<types::global_dof_index> block_sizes = {nqf, np};

    //DynamicSparsityPattern dynamic_sparsity_pattern(dof_handler.n_dofs(), dof_handler.n_dofs());
    BlockDynamicSparsityPattern dynamic_sparsity_pattern(block_sizes, block_sizes);
    DoFTools::make_sparsity_pattern(dof_handler, dynamic_sparsity_pattern);
    constraints.condense(dynamic_sparsity_pattern);

    sparsity_pattern.copy_from(dynamic_sparsity_pattern);
    system_matrix.reinit(sparsity_pattern);

    solution.reinit(block_sizes);
    system_rhs.reinit(block_sizes);

    std::cout<<"setup_system() completed"<<std::endl;
}

// Assemble the system
void Darcyeq::assemble_system() {

    QGauss<2> quadrature_formula(fe.degree + 1);
    QGauss<1> quadrature_formula_face(fe.degree + 1);

    FEValues<2> fe_values (fe, quadrature_formula, update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<2> fe_values_face (fe, quadrature_formula_face, update_values | update_normal_vectors | update_quadrature_points | update_JxW_values);

    // Define Trapezoidal-midpoint quadrature
    // x direction
    std::vector<Point<2>> quad_points_x(2);
    std::vector<double> quad_weights_x(2);
    quad_points_x[0][0] = 0.0; quad_points_x[0][1] = 0.5; 
    quad_points_x[1][0] = 1.0; quad_points_x[1][1] = 0.5;
    quad_weights_x[0] = 0.5; quad_weights_x[1] = 0.5;
    
    // y direction
    std::vector<Point<2>> quad_points_y(2);
    std::vector<double> quad_weights_y(2);
    quad_points_y[0][0] = 0.5; quad_points_y[0][1] = 0.0; 
    quad_points_y[1][0] = 0.5; quad_points_y[1][1] = 1.0; 
    quad_weights_y[0] = 0.5; quad_weights_y[1] = 0.5;
    
    Quadrature<2> quadrature_formula_TMx (quad_points_x, quad_weights_x);
    FEValues<2> fe_values_TMx (fe, quadrature_formula_TMx, update_values | update_gradients | update_quadrature_points | update_JxW_values);
    
    Quadrature<2> quadrature_formula_TMy (quad_points_y, quad_weights_y);
    FEValues<2> fe_values_TMy (fe, quadrature_formula_TMy, update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // Cell matrix and cell rhs vector
    FullMatrix<double> cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
    Vector<double> cell_rhs(fe.dofs_per_cell);

    FullMatrix<double> preconditioner_cell_matrix(fe.dofs_per_cell, fe.dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(fe.dofs_per_cell);

    // Permeability and rhs function objects
    K_inverse K_inverse_object;
    rhs_f rhs_f_object;
    pressure_boundary pressure_boundary_object;

    std::vector<Tensor<2,2>> K_inverse_values(quadrature_formula.size());
    std::vector<Tensor<2,2>> K_inverse_values_TMx (quadrature_formula_TMx.size());
    std::vector<Tensor<2,2>> K_inverse_values_TMy (quadrature_formula_TMy.size());

    double rhs_f_values; 
    std::vector<double> pressure_boundary_values(quadrature_formula_face.size());

    // To store previous iterate values
    std::vector<double> previous_iterate_values (quadrature_formula.size());

    const FEValuesExtractors::Scalar p(2); // dimension = 2
    const FEValuesExtractors::Vector qf(0);
    const FEValuesExtractors::Scalar qfx(0);
    const FEValuesExtractors::Scalar qfy(1);

    // Loop over cells
    for(const auto &cell : dof_handler.active_cell_iterators()) {

        fe_values.reinit(cell);
        cell_matrix = 0.0;
        cell_rhs = 0.0;

        // Get permeability values
        K_inverse_object.value_list(fe_values.get_quadrature_points(), K_inverse_values);

        for (unsigned int q = 0; q < quadrature_formula.size(); ++q) { 

            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) { 
                
                Tensor<1,2> psi_i = fe_values[qf].value(i, q);
                double psi_div_i = fe_values[qf].divergence(i, q);
                double eta_i = fe_values[p].value(i, q);

                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {

                    Tensor<1,2> psi_j = fe_values[qf].value(j, q);
                    double psi_div_j = fe_values[qf].divergence(j, q);
                    double eta_j = fe_values[p].value(j, q);

                    cell_matrix(i,j) += K_inverse_values[q] * psi_i * psi_j * fe_values.JxW(q);
                    cell_matrix(i,j) += -(eta_j * psi_div_i) * fe_values.JxW(q);
                    cell_matrix(i,j) += -(psi_div_j * eta_i) * fe_values.JxW(q);
                }
            }

            // External function values
            rhs_f_values = rhs_f_object.value(fe_values.quadrature_point(q));

            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {

                double eta_i = fe_values[p].value(i, q);

                cell_rhs(i) += -(rhs_f_values * eta_i) * fe_values.JxW(q);
            }
        }

        // Trapezoidal rule implementation 
        /*
        fe_values_TMx.reinit(cell);
        fe_values_TMy.reinit(cell);

        K_inverse_object.value_list(fe_values_TMx.get_quadrature_points(), K_inverse_values_TMx);
        K_inverse_object.value_list(fe_values_TMy.get_quadrature_points(), K_inverse_values_TMy);

        // Trapezoidal midpoint rule
        for (unsigned int q = 0; q < quadrature_formula_TMx.size(); ++q) {
            
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) { 
                Tensor<1,2> psi_x_i = fe_values_TMx[qf].value(i, q);

                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {
                    Tensor<1,2> psi_x_j = fe_values_TMx[qf].value(j, q);

                    if (std::abs(psi_x_i[1] - 0.0) < 1.e-12) {
                        cell_matrix(i,j) += (K_inverse_values_TMx[q] * psi_x_i * psi_x_j) * fe_values_TMx.JxW(q);
                    }
                }
            }
        }

        for (unsigned int q = 0; q < quadrature_formula_TMy.size(); ++q) {
            
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) { 
                Tensor<1,2> psi_y_i = fe_values_TMy[qf].value(i, q);

                for (unsigned int j = 0; j < fe.dofs_per_cell; ++j) {
                    Tensor<1,2> psi_y_j = fe_values_TMy[qf].value(j, q);

                    if (std::abs(psi_y_i[0] - 0.0) < 1.e-12) {
                        cell_matrix(i,j) += (K_inverse_values_TMy[q] * psi_y_i * psi_y_j) * fe_values_TMy.JxW(q);
                    }
                }
            }
        }
       */ 
    
        // Boundary term pressure
        for (unsigned int face_n = 0; face_n < GeometryInfo<2>::faces_per_cell; ++face_n) {
            if ( ( ( cell -> face(face_n) -> at_boundary() ) && ( cell -> face(face_n) -> boundary_id() == 3 ) ) || ( ( cell -> face(face_n) -> at_boundary() ) && ( cell -> face(face_n) -> boundary_id() == 1 ) ) ) {
        
                fe_values_face.reinit(cell, face_n);
                pressure_boundary_object.value_list(fe_values_face.get_quadrature_points(), pressure_boundary_values);
            
                for (unsigned int q = 0; q < quadrature_formula_face.size(); ++q) {
            
                    Tensor<1,2> normal = fe_values_face.normal_vector(q);
                
                    for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) {
                    
                        Tensor<1,2> psi_i_face = fe_values_face[qf].value(i,q);
                    
                        cell_rhs(i) += -(pressure_boundary_values[q] * psi_i_face * normal) * fe_values_face.JxW(q); // Dirichlet pressure conditions
                    }
                }
            }
        }

    (*cell).get_dof_indices (local_dof_indices);
    constraints.distribute_local_to_global(cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
    } // Cell loop

    std::cout<<"assemble_system() completed"<<std::endl;
}

// Solve the system
void Darcyeq::solve() {

    std::cout<<"Solving system using SparseDirectUMFPACK solver"<<std::endl;

    SparseDirectUMFPACK A_direct;
    A_direct.initialize(system_matrix);

    A_direct.vmult(solution, system_rhs);

    constraints.distribute(solution);
    
}

// Output results
void Darcyeq::output_results() {
    
    std::cout<<"Outputting results"<<std::endl;

    std::vector<std::string> solution_names;

    solution_names.push_back("Flux_x");
    solution_names.push_back("Flux_y");
    solution_names.push_back("Pressure");

    DataOut<2> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(dof_handler, solution, solution_names);
    data_out.build_patches();

    const std::string filename = "output/solution" + Utilities::int_to_string (0, 3) + ".vtu";

    std::ofstream output(filename);
    data_out.write_vtu(output);

    // Plot as Vectors
    std::vector<std::string> solution_names_vec(2,"Flux");
    solution_names_vec.emplace_back("Pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation (2, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    data_out.clear_data_vectors();
    data_out.add_data_vector(solution, solution_names_vec, DataOut<2>::type_dof_data, data_component_interpretation);
    data_out.build_patches();

    const std::string filename_vec = "output/solution_vec" + Utilities::int_to_string (0, 3) + ".vtu";

    std::ofstream output_vec(filename_vec);
    data_out.write_vtu(output_vec);

    std::cout<<"Results saved to: "<<filename<<std::endl;
    std::cout<<"output_results() completed"<<std::endl;
}

// Run
void Darcyeq::run() {
    make_grid();
    setup_system();
    assemble_system();
    solve();
    output_results();
}

// Main function
int main() {

    Timer timer;
    std::cout<<"Simulation started"<<std::endl;

    try {
        Darcyeq Darcyeq_obj;
        Darcyeq_obj.run();
    }

    catch(std::exception &exc) {
        std::cerr<<std::endl;
        std::cerr<<"Exception on processing!"<<std::endl
            <<exc.what()<<std::endl    
            <<"Aborting!"<<std::endl;

        return 1;
    }

    catch(...) {
        std::cerr<<"Unknown exception!"<<std::endl
            <<"Aborting!"<<std::endl;
        
        
        return 1;
    }

    timer.stop();

    std::cout<<"Time taken: ("<<timer.cpu_time()<<"s)"<<std::endl;
    std::cout<<"* ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- *"<<std::endl;

    return 0;
}