#ifndef SPC_U_SPC_H_
#define SPC_U_SPC_H_

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <tuple>

#include "macros.h"
#include "u_label.h"

namespace spc {
class USPC {
    public:

        USPC() = default;
        USPC(const USPC&) = delete;
        USPC& operator=(const USPC&) = delete;

    protected:
        void OrderRank() {
            for (uint32_t i = 0; i < n_; ++i) {
                rank_[order_[i]] = i;
            }
        }

        uint32_t n_;
        Graph G_;
        Label dL_, cL_;
    
        std::vector<uint32_t> order_;
        std::vector<uint32_t> rank_;

};

class USPCBasic final: private USPC {
    public:

    USPCBasic() = default;
    USPCBasic(const USPCBasic&) = delete;
    USPCBasic& operator=(const USPCBasic&) = delete;

    std::pair<uint32_t, uint64_t> BFS_SPC(const Graph& const_graph, uint32_t s, uint32_t t);

};

class USPCIndex final: private USPC {
    public:
        enum class OrderScheme {
            kDegree,
            kInvalid
        };

        USPCIndex() = default;
        USPCIndex(const USPCIndex&) = delete;
        USPCIndex& operator=(const USPCIndex&) = delete;

        void BuildIndex(const Graph& const_graph);
        void IndexMerge();
        uint64_t IndexWrite(const std::string& filename);

        void set_os(const OrderScheme os) { os_ = os; }

    private:
        uint32_t Distance(const std::vector<uint32_t>& dLu,
                        const std::vector<LabelEntry>& dLv) const;


        void DegreeOrder(const Graph& graph);
        void InvalidOrder(const Graph&) {
            ASSERT_INFO(false, "invalid ordering");
        }

        std::map<OrderScheme, void (USPCIndex::*)(const Graph&)> of_ = {
            {OrderScheme::kDegree,  &USPCIndex::DegreeOrder},
            {OrderScheme::kInvalid, &USPCIndex::InvalidOrder}
        };

        OrderScheme os_ = OrderScheme::kInvalid;
};

class USPCQuery final: private USPC {
    public:
        USPCQuery() = default;
        USPCQuery(const USPCQuery&) = delete;
        USPCQuery& operator=(const USPCQuery&) = delete;

        std::pair<uint32_t, uint64_t> Count(uint32_t v1, uint32_t v2) const;
        std::pair<uint32_t, uint64_t> bi_BFS_Count(Graph& graph, uint32_t v1, uint32_t v2);

        void IndexRead(const std::string& filename);
        void IndexRead_UPD(const std::string& filename);
        void UpdateGraph(Graph& graph, uint32_t v1, uint32_t v2, char upd_type);
        void print_Label(uint32_t v);

};

class USPCUpdate final: private USPC {
    public:

        USPCUpdate() = default;
        USPCUpdate(const USPCUpdate&) = delete;
        USPCUpdate& operator=(const USPCUpdate&) = delete;

        void IndexRead(const std::string& filename);
        uint64_t IndexWrite(const std::string& filename);

        std::pair<uint32_t, uint64_t> Count(uint32_t v1, uint32_t v2) const;

        std::tuple<uint32_t, size_t, uint32_t, size_t, uint32_t, uint32_t, uint32_t> Inc_SPC(uint32_t a, uint32_t b);
        std::tuple<uint32_t, uint32_t, uint32_t> Inc_BFS(uint32_t hub, uint32_t ab, uint32_t d, uint64_t c);

        std::tuple<uint32_t,uint32_t,size_t,size_t,size_t,size_t,uint32_t,uint32_t,uint32_t,uint32_t> Dec_SPC(uint32_t a, uint32_t b);
        std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> Update_hub(uint32_t hub, 
            const std::vector<int>& Aff_list, const std::vector<uint32_t>& Affs,
            const std::vector<uint32_t>& Recs, int is_hub);
        std::tuple<int, uint32_t, uint32_t, uint32_t, uint32_t> Fast_update(uint32_t a, uint32_t b, 
            const std::vector<int>& Aff_list, const std::vector<uint32_t>& AffA, const std::vector<uint32_t>& AffB,
            const std::vector<uint32_t>& RecA, std::vector<uint32_t> RecB);
        std::tuple<uint32_t, uint64_t, uint32_t, uint64_t, uint32_t> Query_Search(uint32_t h, uint32_t v);
        uint32_t Query_Distance(uint32_t hub, uint32_t v);

    private:
        std::pair<uint32_t, size_t> Distance(const std::vector<uint32_t>& dLu,
                        const std::vector<LabelEntry>& dLv, uint32_t hub) const;

        std::pair<uint32_t, uint32_t> FastDistanceCount(const std::vector<std::pair<uint32_t, uint32_t>>& dLu,
                        const std::vector<LabelEntry>& dLv) const;

        uint32_t FastDistance(const std::vector<uint32_t>& dLu,
                        const std::vector<LabelEntry>& dLv) const;

        std::pair<uint32_t, uint64_t> BFS_SPC(const Graph& const_graph, uint32_t s, uint32_t t);

};

}

#endif
