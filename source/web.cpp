//  This file is part of Artificial Ecology for Chemical Ecology Project
//  Copyright (C) Emily Dolson, 2021.
//  Released under MIT license; see LICENSE

#include <iostream>

#include "emp/prefab/ConfigPanel.hpp"
#include "emp/web/TextArea.hpp"
#include "emp/web/web.hpp"
#include "emp/web/d3/d3_init.hpp"
#include "emp/web/d3/selection.hpp"
#include "emp/web/d3/scales.hpp"
#include "emp/web/d3/svg_shapes.hpp"
#include "emp/web/js_utils.hpp"
#include "emp/datastructs/Graph.hpp"

#include "chemical-ecology/config_setup.hpp"
#include "chemical-ecology/a-eco.hpp"
#include "chemical-ecology/ExampleConfig.hpp"

namespace UI = emp::web;

UI::Document doc("emp_base");
UI::Document cfg_doc("config_panel_div");
// UI::Document viz_doc("grid_visualization");
UI::Canvas canvas;
// UI::Canvas viz_canvas;

emp::vector<D3::LinearScale> square_colors;
emp::vector<std::string> colors;
D3::LinearScale edge_color;
D3::SequentialScale fitness_scale;
D3::SequentialPowScale heredity_scale;


chemical_ecology::Config cfg;
AEcoWorld world;
emp::WeightedGraph interactions;
D3::Selection heatmap;

double max_fitness = 500000;

void ResetScales() {
  EM_ASM({
    emp_i.__outgoing_array = d3.quantize(d3.interpolateHcl("#60c96e", "#4d4193"), $0);
  }, cfg.N_TYPES());
  colors.resize(0);
  emp::pass_vector_to_cpp(colors);
  std::cout << emp::to_string(colors) << std::endl;

  square_colors.resize(cfg.N_TYPES());
  for (int i=0; i < cfg.N_TYPES(); i++){
    square_colors[i].SetDomain(emp::array<double, 2>({0, (double)cfg.MAX_POP()}));
    square_colors[i].SetInterpolate("interpolateLab");
    square_colors[i].SetRange(emp::array<std::string, 2>({"white", colors[i]}));
  }

  edge_color.SetDomain(emp::array<double, 3>({cfg.INTERACTION_MAGNITUDE() * -1, 0, cfg.INTERACTION_MAGNITUDE()}));
  edge_color.SetInterpolate("interpolateLab");
  edge_color.SetRange(emp::array<std::string, 3>({"red", "grey", "blue"}));

  fitness_scale.SetDomain(emp::array<double, 2>({0, max_fitness}));
  fitness_scale.SetInterpolator("interpolateViridis");
  // fitness_scale.SetRange(emp::array<std::string, 2>({"red", "grey", "blue"}));
  heredity_scale.SetDomain(emp::array<double, 2>({0.0, 1.0}));
  heredity_scale.SetInterpolator("interpolateViridis");
  heredity_scale.SetExponent(1.5);
  
}

void DrawWorldCanvas() {
  // UI::Canvas canvas = doc.Canvas("world_canvas");
  canvas.Clear("gray");

  const size_t world_x = cfg.WORLD_X();
  const size_t world_y = cfg.WORLD_Y();
  const double canvas_x = (double) canvas.GetWidth();
  const double canvas_y = (double) canvas.GetHeight();
  const double n_type_sqrt = sqrt(cfg.N_TYPES());
  // std::cout << "x: " << world_x << " y: " << world_y << std::endl;
  const double org_x = canvas_x / (double) world_x;
  const double org_y = canvas_y / (double) world_y;
  const double org_r = emp::Min(org_x, org_y) / 2.0;
  const double type_x = org_x/n_type_sqrt;
  const double type_y = org_y/n_type_sqrt;

  for (size_t y = 0; y < world_y; y++) {
    for (size_t x = 0; x < world_x; x++) {
      const size_t org_id = y * world_x + x;
      if (org_id >= cfg.WORLD_X() * cfg.WORLD_Y()) {
        break;
      }
      const size_t cur_x = org_x * (0.5 + (double) x);
      const size_t cur_y = org_y * (0.5 + (double) y);
      canvas.Rect(x*org_x, y*org_y, org_x, org_y, "#FFFFFF", "black", 10);
      
      for (size_t type = 0; type < cfg.N_TYPES(); type++) {
        int type_x_id = type % (int)n_type_sqrt;
        int type_y_id = type / (int)n_type_sqrt;
        double x_pos = x*org_x + type_x_id*type_x;
        double y_pos = y*org_y + type_y_id*type_y;
        int count = world.world[org_id][type];
        canvas.Rect(x_pos, y_pos, type_x, type_y, square_colors[type].ApplyScale<std::string>(count), "white", .01);        
      }
     }
  }
}

// void DrawWorldViz() {
//   // UI::Canvas canvas = doc.Canvas("world_canvas");
//   viz_canvas.Clear("gray");

//   const size_t world_x = cfg.WORLD_X();
//   const size_t world_y = cfg.WORLD_Y();
//   const double viz_canvas_x = (double) viz_canvas.GetWidth();
//   const double viz_canvas_y = (double) viz_canvas.GetHeight();
//   const double n_type_sqrt = sqrt(cfg.N_TYPES());

//   const double org_x = viz_canvas_x / (double) world_x;
//   const double org_y = viz_canvas_y / (double) world_y;
//   const double org_r = emp::Min(org_x, org_y) / 2.0;
//   const double type_x = org_x/n_type_sqrt;
//   const double type_y = org_y/n_type_sqrt;

//   for (size_t y = 0; y < world_y; y++) {
//     for (size_t x = 0; x < world_x; x++) {
//       const size_t org_id = y * world_x + x;
//       if (org_id >= cfg.WORLD_X() * cfg.WORLD_Y()) {
//         break;
//       }
//       const size_t cur_x = org_x * (0.5 + (double) x);
//       const size_t cur_y = org_y * (0.5 + (double) y);
//       CellData data = world.GetFitness(org_id);
//       if (data.equilib_growth_rate > max_fitness) {
//         max_fitness = data.equilib_growth_rate;
//         fitness_scale.SetDomain(emp::array<double, 2>({0, max_fitness}));
//       }
//       // std::cout << data.heredity << std::endl;
//       viz_canvas.Rect(x*org_x, y*org_y, org_x/2, org_y, fitness_scale.ApplyScale<std::string>(data.equilib_growth_rate), "black", 5);
//       viz_canvas.Rect(x*org_x + org_x/2, y*org_y, org_x/2, org_y, heredity_scale.ApplyScale<std::string>(data.heredity), "black", 5);      

      
//       // std::cout << org_id << ": " << data.fitness << std::endl;
//       // viz_canvas.CenterText({x*org_x + org_x/2, y*org_y + org_y/2}, emp::to_string(data.equilib_growth_rate) + "\n" + emp::to_string((int)data.heredity), "black");
//      }
//   }


  // Add a plus sign in the middle.
  // const double mid_x = org_x * world_x / 2.0;
  // const double mid_y = org_y * world_y / 2.0;
  // const double plus_bar = org_r * world_x;
  // canvas.Line(mid_x, mid_y-plus_bar, mid_x, mid_y+plus_bar, "#8888FF");
  // canvas.Line(mid_x-plus_bar, mid_y, mid_x+plus_bar, mid_y, "#8888FF");

  // doc.Text("ud_text").Redraw();
// }

EM_JS(void, LinkArc, (), {
  d = emp_i.__incoming_array;
  const r = Math.hypot(d[1][0] - d[0][0], d[1][1] - d[0][1]);
  emp.PassStringToCpp(`
    M${d[0][0]},${d[0][1]}
    A${r},${r} 0 0,1 ${d[1][0]},${d[1][1]}
  `);
});

std::string GetLinkArc(emp::array<emp::array<double, 2>, 2> & d) {
  emp::pass_array_to_javascript(d);
  LinkArc();
  return emp::pass_str_to_cpp();
}


void DrawGraph(emp::WeightedGraph g, std::string canvas_id, double radius = 150) {
    D3::Selection s = D3::Select(canvas_id);
    D3::Selection defs = D3::Select("#arrow_defs");

    std::function<std::string(std::string)> MakeArrow = [&defs](std::string color) {
        std::string id = "arrow_" + color;
        emp::remove_chars(id, "#,() ");

        defs.Append("svg:marker")
            .SetAttr("id", id)
            .SetAttr("viewBox", "0 -5 10 10")
            .SetAttr("refX", 5) // This sets how far back it sits, kinda
            .SetAttr("refY", 0)
            .SetAttr("markerWidth", 9)
            .SetAttr("markerHeight", 9)
            .SetAttr("orient", "auto")
            .SetAttr("markerUnits", "userSpaceOnUse")
            .Append("svg:path")
            .SetAttr("d", "M0,-5L10,0L0,5")
            .SetStyle("fill", color);
        
        return "url(#" + id + ")";
    };

    double theta = 0;
    double inc = 2*3.14 / g.GetNodes().size();
    s.SelectAll("circle").Remove();
    for (emp::Graph::Node n : g.GetNodes()) {
        double cx = radius + cos(theta) * radius;
        double cy = radius + sin(theta) * radius;
        theta += inc;
        D3::Selection new_circle = s.Append("circle");
        new_circle.SetAttr("r", 10)
                  .SetAttr("cx", cx)
                  .SetAttr("cy", cy);
                  // .BindToolTipMouseover(node_tool_tip);
        // std::cout << n.GetLabel() << std::endl;
        EM_ASM_ARGS({emp_d3.objects[$0].datum(UTF8ToString($1));}, new_circle.GetID(), n.GetLabel().c_str());
    }

    D3::LineGenerator l;
    s.SelectAll("path").Remove();
    auto weights = g.GetWeights();
    for (int i = 0; i < weights.size(); ++i) {
        for (int j = 0; j < weights[i].size(); ++j) {
            // std::cout << weights[i][j] << std::endl;
            if (abs(weights[i][j]) > .0001) {
                double cxi = radius + cos(i*inc) * radius;
                double cxj = radius + cos(j*inc) * radius;
                double cyi = radius + sin(i*inc) * radius;
                double cyj = radius + sin(j*inc) * radius;
                std::string color = edge_color.ApplyScale<std::string>(weights[i][j]);
                // std::cout << color << std::endl;
                emp::array<emp::array<double, 2>, 2> data({emp::array<double, 2>({cxj, cyj}), emp::array<double, 2>({cxi, cyi})});
                D3::Selection n = s.Append("path");
                // n.SetAttr("d", l.Generate(data))
                n.SetAttr("d", GetLinkArc(data))
                 .SetStyle("stroke", color)
                 .SetStyle("fill", "none")
                //  .SetStyle("stroke-width", edge_width.ApplyScale(weights[i][j]))
                 .SetAttr("marker-end", MakeArrow(color));
                //  .BindToolTipMouseover(edge_tool_tip);

                EM_ASM_ARGS({emp_d3.objects[$0].datum($1);}, n.GetID(), weights[i][j]);
                //  std::cout << MakeArrow(color) << std::endl;
                // s.Append("path").SetAttr("d", "M"+emp::to_string(cxi)+","+emp::to_string(cyi)+"L"+emp::to_string(cxj)+","+emp::to_string(cyj)).SetStyle("stroke", color).SetStyle("stroke-width", weights[i][j]*10);
            }
        }
    }

    // s.SetupToolTip(edge_tool_tip);
    // s.SetupToolTip(node_tool_tip);

}

struct InteractionNode {
  EMP_BUILD_INTROSPECTIVE_TUPLE(
    int, x,
    int, y,
    double, w
  )
};


EM_JS(double, GetNewWeight, (double old), {
  var new_weight = prompt("New weight:", old);
  if (new_weight == null || new_weight == "") {
    return old;
  }
  return +new_weight;
});

emp::WeightedGraph MakeGraph() {
  emp::WeightedGraph interaction_graph;

  interaction_graph.Resize(cfg.N_TYPES());
  emp::vector<emp::vector<double> > int_matrix = world.GetInteractions();

  for (int i = 0; i < cfg.N_TYPES(); i++) {
    for (int j = 0; j < cfg.N_TYPES(); j++) {
      interaction_graph.AddEdge(i, j, int_matrix[i][j]);
    }
  }

  return interaction_graph;
}


void DrawInteractionMatrix(emp::WeightedGraph & g, std::string canvas_id, int width = 35) {
    emp::vector<InteractionNode> interaction_vec;
    interaction_vec.resize(cfg.N_TYPES()*cfg.N_TYPES());

    for (int x = 0; x < cfg.N_TYPES(); x++) {
      for (int y = 0; y < cfg.N_TYPES(); y++) {
        InteractionNode & n = interaction_vec[y + x * cfg.N_TYPES()];
        n.x(x);
        n.y(y);
        n.w(g.GetWeight(x,y));
      }
    }
    
    D3::Selection s = D3::Select(canvas_id);
    D3::Selection enter = s.SelectAll("rect")
     .Data(interaction_vec)
     .EnterAppend("rect");

    heatmap = enter.Merge(heatmap);

    heatmap
     .SetAttr("x", [width](InteractionNode n){return n.x()*width;})
     .SetAttr("y", [width](InteractionNode n){return n.y()*width;})
     .SetAttr("width", [width](InteractionNode n){return width;})
     .SetAttr("height", [width](InteractionNode n){return width;})
     .SetStyle("fill", [](InteractionNode n){return edge_color.ApplyScale<std::string>(n.w());})
     .On("click", [&g, canvas_id, width](InteractionNode n, int id){
      //  std::cout << "Clicked" << std::endl;
       double new_weight = GetNewWeight(n.w());
       g.AddEdge(n.x(), n.y(), new_weight);
       world.SetInteraction(n.x(), n.y(), new_weight);
      //  std::cout << g.GetWeight(n.x(), n.y()) << std::endl;
      //  EM_ASM({
      //    var d = .datum();
      //    d.w = $1;
      //    emp_d3.objects[$0].datum(d);
      //    emp_d3.objects[$0].style("fill", emp_d3.objects[$1]($0));
      //  }, id, new_weight, edge_color.GetID());
      DrawInteractionMatrix(g, canvas_id, width);
      interactions = MakeGraph();
      DrawGraph(interactions, "#interaction_network");
     });

    emp::vector<emp::vector<double> > ints = world.GetInteractions();

    for (int i = 0; i < cfg.N_TYPES(); i++) {
      std::cout << emp::join(ints[i], ",") << std::endl;
    }

}

emp::Ptr<UI::Animate> anim;


int main()
{
  // n_objects;
  // doc << "<h1>Hello, browser!</h1>";
  canvas = doc.AddCanvas(800, 800, "world_canvas");
  canvas.SetHeight(800);
  canvas.SetWidth(800);
  // doc << canvas;

  // viz_canvas = viz_doc.AddCanvas(600, 600, "world_viz");
  // viz_canvas.SetHeight(600);
  // viz_canvas.SetWidth(600);


  anim.New([](){ world.Update(world.GetTime()); DrawWorldCanvas(); /*DrawWorldViz();*/ }, canvas );

  // Set up a configuration panel for web application
  setup_config_web(cfg);
  cfg.Write(std::cout);
  emp::prefab::ConfigPanel example_config_panel(cfg);
  // example_config_panel.ExcludeSetting("SEED");
  // example_config_panel.SetRange("SEED", "-1", "100", "1");
  cfg_doc << example_config_panel;
  cfg_doc << anim->GetToggleButton("but_toggle");
  // cfg_doc << UI::Button([](){DrawWorldCanvas(true);}, "Calculate fitness");

  // An example to show how the Config Panel could be used
  // to control the color of text in an HTML text area
  // UI::TextArea example_area("example_area");
  // example_area.SetSize(cfg.SIZE(), cfg.SIZE());
  // example_config_panel.SetOnChangeFun([](const std::string & setting, const std::string & value){

  // });

  world.Setup(cfg);
  // world.Run();
  ResetScales();
  DrawWorldCanvas();
  // DrawWorldViz();
  interactions = MakeGraph();
  DrawGraph(interactions, "#interaction_network");
  DrawInteractionMatrix(interactions, "#interaction_matrix");
}
