import React from "react";

import { List } from "@mui/material";

const SearchResultList = (props) => {
  const ListItemWrapper = props.listItemWrapper;

  return (
    <List style={{ display: "flex", flexDirection: "column" }}>
      {props.searchedData.map((searchedItem) => (
        <ListItemWrapper key={searchedItem.id} item={searchedItem} />
      ))}
    </List>
  );
};

export default SearchResultList;
