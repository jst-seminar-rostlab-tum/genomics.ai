import { ExpandLess, ExpandMore } from '@mui/icons-material';
import {
  Box, Collapse, IconButton, Tooltip, Typography,
} from '@mui/material';
import React, { useState } from 'react';

/**
 * A sidepanel for the GeneMapper to show additional information
 * that can be collaped to either the right or left side.
 * @param title Heading to show at the top of the sidepanel
 * @param collapseToRight Set to true if panel should be collapsing to the right
 * @param children Child nodes inside the panel
 */
function Sidepanel({ title, collapseToRight, children }) {
  const [open, setOpen] = useState(true);

  const expandIconAngle = collapseToRight ? 90 : -90;
  const toolTipPlacement = collapseToRight ? 'left' : 'right';
  const flexDirection = collapseToRight ? 'row-reverse' : 'row';

  return (
    <Box sx={{
      display: 'flex',
      flexDirection: 'column',
      alignItems: 'strech',
      height: '100%',
    }}
    >
      {open
        ? (
          <Box sx={{ display: 'flex', alignItems: 'center', flexDirection }}>
            <Typography variant="h6" style={{ flexGrow: 1 }}>{title}</Typography>
            <IconButton onClick={() => setOpen(false)}>
              <ExpandLess fontSize="large" sx={{ transform: `rotate(${expandIconAngle}deg)` }} />
            </IconButton>
          </Box>
        )
        : (
          <Tooltip title={title} placement={toolTipPlacement}>
            <IconButton onClick={() => setOpen(true)}>
              <ExpandMore fontSize="large" sx={{ transform: `rotate(${expandIconAngle}deg)` }} />
            </IconButton>
          </Tooltip>
        )}
      <Collapse in={open} orientation="horizontal" sx={{ overflowY: 'auto', overflowX: 'hidden', pr: 2 }}>
        {children}
      </Collapse>
    </Box>
  );
}

export default Sidepanel;
