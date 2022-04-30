import React from "react";

import { Button, Chip } from "@mui/material";

import Avatars from "components/Avatars";
import SearchCard from "../SearchCard";
import LabeledLink from "../LabeledLink";

import Tag from "components/Tag";
import CustomButton from "components/CustomButton"

// Card to display search result for a single team
const TeamCard = ({ item: team }) => {
  return (
    <SearchCard
      // action={<Button variant="contained">Join</Button>}
      action={<CustomButton type="primary">Join</CustomButton>}
      title={team.name}
      primary={
        //  <Tag content={team.visibility} variant="primary-default" />
        <Chip label={team.visibility} color="primary" size="small" />
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
          {team.institution && (
            <LabeledLink label={"Institution"} content={team.institution} />
          )}
        </React.Fragment>
      }
    />
  );
};

export default TeamCard;
