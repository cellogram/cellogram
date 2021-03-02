////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <zebrafish/Logger.hpp>

////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

using zebrafish::logger;

namespace {

struct PropertyEditorItem {

static void AppendMarkerRecordItem(
    const char *prefix, 
    const int uid, 
          Eigen::MatrixXd &marker_4D, 
    const zebrafish::OptimDepthInfo_t &depthInfo, 
    const zebrafish::DepthSearchFlag_t &flag) {

    std::string postfix = "";
    if (flag == zebrafish::DepthSearchFlag_t::Success)
        postfix = "";
    else if (flag == zebrafish::DepthSearchFlag_t::InvalidEnergy)
        postfix = " [Ener]";
    else if (flag == zebrafish::DepthSearchFlag_t::SecondDerivative)
        postfix = " [Deri]";
    else
        postfix = " [Unkn]";

    ImGui::PushID(uid);
    ImGui::AlignTextToFramePadding();
    bool nodeOpen = ImGui::TreeNode("Object", "%s %u%s", prefix, uid, postfix.c_str());
    ImGui::NextColumn();
    ImGui::AlignTextToFramePadding();

    // no text here

    ImGui::NextColumn();

    if (nodeOpen) {
        static const std::vector<std::string> itemName{"X (col)", "Y (row)", "Z", "R", "flag"};
        for (int i = 0; i < 4; i++) {
            ImGui::PushID(i);
            ImGui::AlignTextToFramePadding();
            ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
            ImGui::NextColumn();

            // ImGui::Text("%.3f", markerRecord.loc(uid, i));
            if (ImGui::InputDouble("", &marker_4D(uid, i), 0.0, 0.0, "%.3f")) {

            }

            ImGui::NextColumn();
            ImGui::PopID();
        }

        ImGui::Separator(); ///////////////////////

        ImGui::PushID(4);
        ImGui::AlignTextToFramePadding();
        ImGui::TreeNodeEx(itemName[4].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
        ImGui::NextColumn();
        if (flag == zebrafish::DepthSearchFlag_t::Success)
            ImGui::Text("Success");
        else if (flag == zebrafish::DepthSearchFlag_t::InvalidEnergy)
            ImGui::Text("Invalid Energy");
        else if (flag == zebrafish::DepthSearchFlag_t::SecondDerivative)
            ImGui::Text("Second Derivative");
        else
            ImGui::Text("Unknown");
        ImGui::NextColumn();
        ImGui::PopID();

        ImGui::Separator(); ///////////////////////
        ImGui::TreePop();
    }

    ImGui::PopID();
}
};  // struct PropertyEditorItem

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// Special module about property editor

void UIState::draw_editor_window() {

    static bool firstTimeReachHere = true;
    if (UIsize.resize || firstTimeReachHere) {
		ImGui::SetNextWindowPos(ImVec2(UIsize.windowWidth-UIsize.rightWidth, UIsize.mainMenuHeight));
		ImGui::SetNextWindowSize(ImVec2(UIsize.rightWidth, UIsize.windowHeight-UIsize.mainMenuHeight-UIsize.imageViewerHeight));
        firstTimeReachHere = false;  // because this is not open by default
	}

	if (!ImGui::Begin("Property Editor", &show_editor_menu)) {
		ImGui::End();
		return;
	}

    ImGui::PushItemWidth(UIsize.rightWidth / 2.0);
    std::vector<std::string> typeName{"Marker"};
    static int propertyListType = 0;
    ImGui::Combo("Property List Type", &propertyListType, typeName);
    ImGui::PopItemWidth();

    ImGui::Separator();
    ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

    switch (propertyListType) {
    case 0:
        // markers (finalized clusters)
        if (state.mesh.C_depth_info_vec.empty()) {
            ImGui::Text("3D marker list is empty");
        } else {
            ImGui::PushStyleVar(ImGuiStyleVar_FramePadding, ImVec2(2, 2));
            ImGui::Columns(2);

            static int maxNumItemDisplayed = 200;
            const int ttlItem = state.mesh.C_depth_info_vec.size();
            const int numItemToDisplay = std::min(maxNumItemDisplayed, ttlItem);
            for (int i = 0; i < numItemToDisplay; i++) {
                PropertyEditorItem::AppendMarkerRecordItem(
                    "Marker", 
                    i, 
                    state.mesh.marker_4D,
                    state.mesh.C_depth_info_vec[i], 
                    state.mesh.dsFlag[i]);
            }

            ImGui::Columns(1);
            if (ttlItem >= maxNumItemDisplayed) {
                ImGui::Text("Only the first %d items will be displayed", maxNumItemDisplayed);
                if (ImGui::Button("Show more")) {
                    maxNumItemDisplayed += 100;
                }
            }
            ImGui::PopStyleVar();
        }
        break;

    default:
        assert(false);
        break;
    }

    ImGui::EndChild();
    ImGui::End();
}

}  // namespace cellogram
