import React from "react";

import { Divider, List } from "@mui/material";
import TeamCard from "./TeamCard";

const SearchResultList = (props) => {
  return (
    <List style={{ display: "flex", flexDirection: "column" }}>
      {props.searchedData.map((searchedItem) => (
        <TeamCard team={searchedItem} />
      ))}
    </List>
  );
};

export default SearchResultList;
