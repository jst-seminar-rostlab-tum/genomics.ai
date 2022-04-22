import React from "react";

import { Button, Chip } from "@mui/material";
import SearchCard from "../SearchCard";
import Avatars from "components/Avatars";

const TeamCard = ({ team }) => {
  return (
    <SearchCard
      id={team.id}
      action={<Button variant="contained">Join</Button>}
      title={team.name}
      primary={
        <React.Fragment>
          <Chip label={team.visibility} color="primary" size="small" />
        </React.Fragment>
      }
      secondary={`updated on ${team.updated}`}
      tertiary={
        <React.Fragment>
          <Avatars
            items={team.members.map(({ name, image }) => {
              return { src: image, alt: name };
            })}
          />
          <Chip
            label={`${team.membersCount} members`}
            variant="outlined"
            size="small"
            sx={{ color: "text.secondary" }}
          />
        </React.Fragment>
      }
    />
  );
};

export default TeamCard;
