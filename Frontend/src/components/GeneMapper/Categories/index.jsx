import {
  Circle, ExpandLess, ExpandMore, FormatColorFill, FormatColorReset, Visibility, VisibilityOff,
} from '@mui/icons-material';
import {
  Box, Collapse, IconButton, ListItem, ListItemButton, Typography,
} from '@mui/material';
import React, { useState } from 'react';

function Category({ title, colored, toggleColored }) {
  const [open, setOpen] = useState(false);

  return (
    <div style={{ display: 'flex', flexDirection: 'column' }}>
      <ListItemButton onClick={() => setOpen(!open)} disableGutters>
        <ListItem
          secondaryAction={(
            <IconButton onClick={(e) => { toggleColored(); e.stopPropagation(); }}>
              {colored ? <FormatColorFill /> : <FormatColorReset />}
            </IconButton>
          )}
        >
          <Typography sx={{ pr: 1 }}>{title}</Typography>
          {open ? <ExpandLess /> : <ExpandMore />}

        </ListItem>
      </ListItemButton>
      <Collapse in={open}>
        <div style={{ display: 'flex', flexDirection: 'column' }}>
          <Value title="short" />
          <Value title="veeeeeeeeeeeeeeeery long" />
          <Value title="medium length" />
        </div>
      </Collapse>
    </div>
  );
}

function Value({ title }) {
  const [visible, setVisible] = useState(true);

  return (
    <Box sx={{ display: 'flex', alignItems: 'center', px: 3 }}>
      <IconButton edge="start" onClick={() => setVisible(!visible)}>
        {visible ? <Visibility /> : <VisibilityOff />}
      </IconButton>
      <Typography sx={{ flexGrow: 1, pr: 3 }}>{title}</Typography>
      <Circle />
    </Box>
  );
}

function GeneMapperCategories({ categories }) {
  const [coloredCategoryIndex, setColoredCategoryIndex] = useState(0);
  const [open, setOpen] = useState(true);

  return (
    <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'start' }}>
      <Box sx={{ display: 'flex', flexDirection: 'column', alignItems: 'strech' }}>
        <Box sx={{ display: 'flex', alignItems: 'center', pr: 1 }}>
          <Typography variant="h6" style={{ flexGrow: 1 }}>Categories</Typography>
          <IconButton size="large">
            <ExpandLess fontSize="large" sx={{ transform: open ? 'rotate(90deg)' : 'rotate(-90deg)' }} />
          </IconButton>
        </Box>
        {categories.map((category, idx) => (
          <Category
            key={category}
            title={category}
            colored={idx === coloredCategoryIndex}
            toggleColored={() => {
              setColoredCategoryIndex(idx);
            }}
          />
        ))}
      </Box>
    </Box>
  );
}

export default GeneMapperCategories;
