import React from "react";

import {
  ListItem,
  Stack,
  Link,
  Typography,
  ListItemAvatar,
  Avatar,
} from "@mui/material";
import { Link as RouterLink } from "react-router-dom";

// Temporary function to generate different colors, change later to titleToColor(...) from feature/institutinOverview branch
function randomColor() {
  let hex = Math.floor(Math.random() * 0xffffff);
  let color = "#" + hex.toString(16);

  return color;
}

// Generic card component to reuse the same structure for the different search card items as Institution, team, ...
const SearchCard = (props) => {
  return (
    <ListItem divider alignItems="flex-start" secondaryAction={props.action}>
      {props.avatar && (
        <ListItemAvatar spacing={10}>
          <Avatar 
            sx={{ bgcolor: randomColor(), width: 45, height: 45 }}
            alt={props.title}
            src={props.avatar}
          />
        </ListItemAvatar>
      )}
      <Stack direction="column" spacing={0.5}>
        {/* primary */}
        <Stack direction="row" spacing={2} alignItems="center">
          <Link
            sx={{ fontSize: "24px" }}
            component={RouterLink}
            to="documentation" // TODO integarte as prop
            underline="hover"
          >
            <b>{props.title}</b>
          </Link>
          {props.primary}
        </Stack>
        {/* secondary */}
        <Typography sx={{ color: "text.secondary" }}>
          {props.secondary}
        </Typography>
        {/* tertiary */}
        <Stack
          direction="row"
          sx={{ color: "text.secondary" }}
          alignItems="flex-end"
          spacing={2}
        >
          {props.tertiary}
        </Stack>
      </Stack>
    </ListItem>
  );
};

export default SearchCard;
