import React from "react";

import { List } from "@mui/material";
import InstitutionCard from "./InstitutionCard";
import TeamCard from "./TeamCard";

const SearchResultList = (props) => {
  const createCard = (item) => {
    switch (props.type) {
      case "Teams":
        return <TeamCard team={item} />;
      case "Institutions":
        return <InstitutionCard institution={item} />;
    }
  };

  return (
    <List style={{ display: "flex", flexDirection: "column" }}>
      {props.searchedData.map((searchedItem) => {
        return <React.Fragment key={searchedItem.id} > {createCard(searchedItem)} </React.Fragment>;
      })}
    </List>
  );
};

export default SearchResultList;
