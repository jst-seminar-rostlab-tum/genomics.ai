import React from "react";

import {
  ListItem,
  List,
  ListItemText,
  Button,
  AvatarGroup,
  Avatar,
  Stack,
  ListItemAvatar,
} from "@mui/material";
import { deepOrange } from "@mui/material/colors";
import SearchCard from "./SearchCard";

// Temporary function to generate diff colors, use later from feature/institutioOverview branch
function randomColor() {
  let hex = Math.floor(Math.random() * 0xffffff);
  let color = "#" + hex.toString(16);

  return color;
}

const SearchResultList = (props) => {
  return (
    <List style={{ display: "flex", flexDirection: "column", maxWidth: 800 }}>
      {props.searchedData.map((searchedItem) => (
        <SearchCard
          key={searchedItem.id}
          action={<Button variant="contained">Join</Button>}
          primary={
            <React.Fragment>
              <div>{searchedItem.name}</div>
              <div>{searchedItem.visibility}</div>
            </React.Fragment>
          }
          secondary={searchedItem.updated}
          tertiary={
            <React.Fragment>
              <AvatarGroup>
                {searchedItem.members.map((member) => (
                  <Avatar
                    sx={{ bgcolor: randomColor(), width: 24, height: 24 }}
                    alt={member.name}
                    src="/static/images/avatar/1.jpg"
                  />
                ))}
              </AvatarGroup>
              <div>{`${searchedItem.membersCount} members`}</div>
            </React.Fragment>
          }
        />

      ))}
    </List>
  );
};

export default SearchResultList;
