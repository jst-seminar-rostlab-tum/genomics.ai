import React from "react";

import { List } from "@mui/material";
import InstitutionCard from "./InstitutionCard";
import TeamCard from "./TeamCard";

const SearchResultListItem = (props) => {
  switch (props.type) {
    case "Teams":
      return <TeamCard team={props.item} />;
    case "Institutions":
      return <InstitutionCard institution={props.item} />;
  }
};

const SearchResultList = (props) => {
  return (
    <List style={{ display: "flex", flexDirection: "column" }}>
      {props.searchedData.map((searchedItem) => {
        return (
          <React.Fragment key={searchedItem.id}>
            <SearchResultListItem type={props.type} item={searchedItem} />
          </React.Fragment>
        );
      })}
    </List>
  );
};

export default SearchResultList;
