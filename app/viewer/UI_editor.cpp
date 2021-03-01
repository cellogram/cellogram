////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <zebrafish/Logger.hpp>

////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

using zebrafish::logger;

namespace {
/*
struct PropertyEditorItem {

    static bool AppendMarkerRecordItem(const char *prefix, int uid, markerRecord_t &markerRecord) {

        ImGui::PushID(uid);
        ImGui::AlignTextToFramePadding();
        bool nodeOpen = ImGui::TreeNode("Object", "%s %u", prefix, uid);
        ImGui::NextColumn();
        ImGui::AlignTextToFramePadding();
        bool res = false;

        // no text here

        ImGui::NextColumn();

        if (nodeOpen) {
            static const std::vector<std::string> itemName{"X (row)", "Y (col)", "Z", "R", "energy", "size"};
            for (int i = 0; i < 4; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();

                // ImGui::Text("%.3f", markerRecord.loc(uid, i));
                if (ImGui::InputDouble("", &markerRecord.loc(uid, i), 0.0, 0.0, "%.3f")) {
                    res = true;
                }

                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////

            for (int i = 4; i <= 5; i++) {
                ImGui::PushID(i);
                ImGui::AlignTextToFramePadding();
                ImGui::TreeNodeEx(itemName[i].c_str(), ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen | ImGuiTreeNodeFlags_Bullet);
                ImGui::NextColumn();
                if (i == 4) {
                    ImGui::Text("%.4f", markerRecord.energy(uid));
                } else {
                    ImGui::Text("%d", markerRecord.size(uid));
                }
                ImGui::NextColumn();
                ImGui::PopID();
            }

            ImGui::Separator(); ///////////////////////
            ImGui::TreePop();
        }

        ImGui::PopID();
        return res; // Whether InputDouble triggered
    }
};
*/
} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// Special module about property editor

} // namespace cellogram
