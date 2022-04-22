import React from "react";
import { AvatarGroup, Avatar } from "@mui/material";
import { grey } from "@mui/material/colors";
import MoreHorizIcon from "@mui/icons-material/MoreHoriz";

// Temporary function to generate different colors, change later to titleToColor(...) from feature/institutinOverview branch
function randomColor() {
  let hex = Math.floor(Math.random() * 0xffffff);
  let color = "#" + hex.toString(16);

  return color;
}

// A component displaying a common group of avatrs
const Avatars = ({
  items,
  width = 24,
  height = 24,
  displayMoreIcon = true,
}) => {
  return (
    <AvatarGroup>
      {items.map((item) => (
        <Avatar
          key={item.src}
          sx={{ bgcolor: randomColor(), width, height }}
          alt={item.alt}
          src={item.src}
        />
      ))}
      { displayMoreIcon && <Avatar sx={{ bgcolor: grey[600], width, height }} src="/dummy.jpg">
        <MoreHorizIcon />
      </Avatar>}
    </AvatarGroup>
  );
};

export default Avatars;
